import sys
import numpy as np

from amino_data import letter_code, structures, Atom
from pdb_util import pdb_line
from structures_codegen import extract_atom_positions


# Maximum distance from any atom in an amino to the peptide nitrogen
MAX_DIST_TO_N = 10.0              # [Å] rounded up from 8.48
ATOM_RADIUS = 0.68                # [Å] atom radius for C, N, O
COLLISION_DIST = 3.*ATOM_RADIUS   # [Å] twice atom radius to form a bond, plus 1.0 safety margin

# Search constants:
P_BACKTRACK = 0.2                 # base probability that we backtrack a given bad chain even farther
SEARCH_WIDTH = 16                 # width of the DFS
PROGRESS_SPACING = 20             # show progress every this many amino acids


# to make sure we can handle proteins with a large number of residues, increase the recursion limit:
sys.setrecursionlimit(2048)


def blks_to_pdb(residues, blks):
    lines = []
    i_atom = 0
    for i_seq, (residue, blk) in enumerate(zip(residues, blks)):
        for atom in blk:
            lines.append(pdb_line(atom, i_atom, 1 + i_seq, residue))
            i_atom += 1
    lines.append("TER " + pdb_line(Atom(" ", [0., 0., 0.], " "), i_atom, len(residues), residue)[4:27])
    return "\n".join(lines)

def blks_to_cg_pdb(residues, blks, output_filename):
    """ Writes a CG PDB file where each particle corresponds to the C-alpha (CA) atom of each residue. """
    lines = []
    i_atom = 1  # Start atom numbering from 1
    for i_seq, (residue, blk) in enumerate(zip(residues, blks)):
        # Find the C-alpha (CA) atom in the block
        ca_atom = next((atom for atom in blk if atom.name == "CA"), None)
        if ca_atom:
            # Write the CA atom as the CG particle
            lines.append(pdb_line(ca_atom, i_atom, 1 + i_seq, residue))
            i_atom += 1
    
    # Add the TER statement at the end of the file
    lines.append(f"TER   {i_atom:>5}")
    
    # Write to the CG PDB file
    with open(output_filename, "w") as output_file:
        output_file.write("\n".join(lines))
        output_file.write("\n")

class PartialStructure:
    """ class representing the partial structure of an unfolded chain.
        keeps track of current displacement and spatial rotation,
        sequence of amino acid choices so far, and allows backtracking
        through the history of building the chain """
    def __init__(self, full_length):
        self.i = 0
        self.len = full_length
        self.cursor = np.zeros((full_length + 1, 3))            # indexed by i
        self.transform = np.zeros((full_length + 1, 3, 3))      # indexed by i
        self.transform[0] = np.identity(3)
        self.aminos = [None for _ in range(full_length)]        # indexed by i - 1
        self.blks = [None for _ in range(full_length)]          # indexed by i - 1
        self.atom_poses = [None for _ in range(full_length)]    # indexed by i - 1
    def _transform_pos(self, pos):
        return self.cursor[self.i] + (self.transform[self.i] @ pos)
    def add_struct(self, amino_struct):
        blk = []
        for atom in amino_struct.atoms:
            transformed_atom = atom._replace(pos = self._transform_pos(atom.pos))
            if atom.name == "N_next": # NOTE: N_next should be the last element in the list!
                self.cursor[self.i + 1] = transformed_atom.pos
            else:
                blk.append(transformed_atom)
        self.transform[self.i + 1] = self.transform[self.i] @ amino_struct.transform
        self.aminos[self.i] = amino_struct
        self.blks[self.i] = blk
        self.atom_poses[self.i] = np.stack([atom.pos for atom in blk]) # np array for faster processing
        self.i += 1
        if self.i == self.len:
            self._terminate_last_amino()
    def _terminate_last_amino(self):
        assert self.i == self.len, "haven't yet reached end of chain"
        last_amino = self.blks[-1]
        (CA_pos, C_pos, O_pos) = extract_atom_positions(["CA", "C", "O"], last_amino)
        OC2_pos = 3*C_pos - (CA_pos + O_pos)
        idx_O = [i for i in range(len(last_amino)) if last_amino[i].name == "O"][0]
        last_amino[idx_O] = last_amino[idx_O]._replace(name = "OC1")
        last_amino.append(Atom("OC2", OC2_pos, "O"))
    def count_collisions(self):
        """ count collisions with the top (most recently added) residue """
        collisions = 0
        i = self.i - 1 # actual index for the top amino
        for j in range(i):
            poses_i, poses_j = self.atom_poses[i], self.atom_poses[j]
            if np.linalg.norm(poses_i[0] - poses_j[0]) > 2*MAX_DIST_TO_N + COLLISION_DIST:
                continue # guaranteed no collisions thanks to triangle inequality
            if j == i - 1: # when checking for collisions with a neighbouring residue, bonded atoms allowed to collide
                poses_i = poses_i[1:] # pull off the N atom
            collisions += (np.linalg.norm(poses_i[:, None] - poses_j[None, :], axis=2) < COLLISION_DIST).sum()
        return collisions
    def rough_radius(self):
        """ compute rough radius based on Nitrogen positions of current aminos """
        poses = self.cursor[:self.i]
        pos_avg = poses.mean(0)
        return np.sqrt(((poses - pos_avg)**2).sum(1).max(0))
    def clone(self):
        """ create a duplicate of this amino acid """
        ans = PartialStructure(self.len)
        for j in range(self.i):
            ans.add_struct(self.aminos[i])
    def pop(self):
        """ roll back by popping off top amino acid """
        self.i -= 1


def pick_random_struct(residue):
    choices = structures[residue]
    return choices[np.random.randint(len(choices))]


def get_smallest_chain(chains):
    largest_radius = 1000000000000
    ans = None
    for j, chain in enumerate(chains):
        r = chain.rough_radius()
        if r < largest_radius:
            largest_radius = r
            ans = j
    return ans


class ProgressBar:
    def __init__(self):
        self.highest_visited = 0
    def register_progress(self, i):
        if i % PROGRESS_SPACING == 0 and i > self.highest_visited:
            sys.stderr.write("reached %d residues\n" % i)
            self.highest_visited = i


def wide_df_search(i, chains, sequence, progress_bar=None):
    """ creates a recursive pile of python generators, it's a little hacky.
        yield value is a list of 'bad chains' that must be readjusted.
        yield value of None means the search has completed. """
    if progress_bar is None:
        progress_bar = ProgressBar()
    if i < len(sequence): # generator recursion!
        sub_search_generator = wide_df_search(i + 1, chains, sequence, progress_bar)
    else:
        def none_generator():
            yield None
            return
        sub_search_generator = none_generator()
    if i > 0: yield list(range(len(chains))) # all sub-chains will initially need updating
    for bad_chains_sub in sub_search_generator:
        if bad_chains_sub is None: # search complete, yay!
            yield None
            return
        while True:
            progress_bar.register_progress(i)
            bad_chains = []
            for j in bad_chains_sub:
                while True:
                    if np.random.rand() < P_BACKTRACK*(i/(i + 1)): # adjust backtrack probability based on how near we are to the start of the sequence
                        if j not in bad_chains:
                            bad_chains.append(j)
                        break
                    chains[j].add_struct(pick_random_struct(sequence[i]))
                    if chains[j].count_collisions() > 0:
                        chains[j].pop() # pop the thing we just added so there's space for our next try
                    else:
                        break # successful addition of amino
            if len(bad_chains) > 0:
                for j in bad_chains:
                    chains[j].pop()
                yield bad_chains # fix these chains on the levels above this one
                bad_chains_sub = bad_chains # they still need to be updated at this level too
            else:
                break


def pdb_chain(sequence, showprogress=True):
    sequence = [letter_code[c] for c in sequence] # convert to 3 letter code
    chains = [PartialStructure(len(sequence)) for _ in range(SEARCH_WIDTH)]
    _ = list(wide_df_search(0, chains, sequence)) # do the search
    sys.stderr.write("generated radii: %s\n" % ",  ".join([str(chain.rough_radius()) for chain in chains]))
    best_chain_idx = get_smallest_chain(chains)
    sys.stderr.write("best: %f\n" % chains[best_chain_idx].rough_radius())
    return blks_to_pdb(sequence, chains[best_chain_idx].blks), chains, best_chain_idx
    #return blks_to_pdb(sequence, chains[best_chain_idx].blks)

if __name__ == "__main__":
    # Get the input sequence and convert it to the 3-letter amino acid code
    raw_sequence = input()  # This will be the input amino acid sequence (e.g., "MASNDYT...")
    sequence = [letter_code[c] for c in raw_sequence]  # Convert to 3-letter code

    # Generate the full atomistic PDB structure and get chains and best chain index
    full_pdb_content, chains, best_chain_idx = pdb_chain(raw_sequence)
    
    # Print or write the full atomistic PDB content
    print(full_pdb_content)
    
    # Now generate the coarse-grained PDB where each particle corresponds to the C-alpha atoms
    cg_pdb_filename = 'cg_fus.pdb'  # Name of the CG PDB file
    blks_to_cg_pdb(sequence, chains[best_chain_idx].blks, cg_pdb_filename)
    print(f"Coarse-grained PDB file saved as {cg_pdb_filename}")

