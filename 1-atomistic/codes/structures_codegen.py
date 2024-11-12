import numpy as np
from collections import defaultdict

from pdb_util import parse_pdb, write_xyz, ResidueExample


# This file is meant for code generation.
# Given a particular pdb file, we pick some representatives of each
# kind of amino acid from the file. We then generate python code
# for a dict that contains these structures.


PDB_PATH = "7sbp.pdb" # can obtain this file here: https://files.rcsb.org/download/7SBP.pdb
REPS_LIST = [ # list of representative residues
    *["A%d" % i for i in range(50, 69)],
    *["A%d" % i for i in range(78, 174)],
    *["A%d" % i for i in range(186, 245)],
    *["A%d" % i for i in range(256, 676)],
    *["A%d" % i for i in range(689, 1000)]]


def tounit(v):
    return v / np.linalg.norm(v)

def get_unit_basis(CA_pos, C_pos):
    """ return unit basis with orientations standardized so that
        CA is on X axis and C is in XY plane
        we assume that N is at the origin.
        CA_pos: (3) - position of CA relative to N
        C_pos:  (3) - position of C relative to N """
    v = tounit(C_pos - CA_pos)
    # make the basis vectors:
    e_0 = tounit(CA_pos)
    e_2 = tounit(np.cross(e_0, v))
    e_1 = tounit(np.cross(e_2, e_0))
    return e_0, e_1, e_2

def extract_atom_positions(atom_name_list, blk):
    ans = [None for name in atom_name_list]
    for atom in blk:
        if atom.name in atom_name_list:
            idx = atom_name_list.index(atom.name)
            assert ans[idx] is None, "duplicate atom in residue block!"
            ans[idx] = atom.pos
    return ans

def normalize_positions(blk):
    """ standardize positions and orientations of atoms in a residue """
    # translate:
    (origin,) = extract_atom_positions(["N"], blk)
    blk = [atom._replace(pos = np.array(atom.pos) - origin) for atom in blk]
    # rotate:
    (CA_pos, C_pos) = extract_atom_positions(["CA", "C"], blk)
    rot_matrix = np.stack(get_unit_basis(CA_pos, C_pos), axis=0)
    return [atom._replace(pos = rot_matrix @ atom.pos) for atom in blk]

def get_transform(blk):
    """ get the transform matrix between successive amino acids.
        ASSUMES: the positions and orientations have already been normalized.
        the result here will be based on the positions of the atoms in the next residue. """
    (N_pos, CA_pos, C_pos) = extract_atom_positions(["N_next", "CA_next", "C_next"], blk)
    # note that we use axis=-1 in this case because we the unit vectors to map to the next basis
    return np.stack(get_unit_basis(CA_pos - N_pos, C_pos - N_pos), axis=-1)


def get_residue_blocks(atoms):
    ans = defaultdict(lambda: [])
    # add regular atoms
    for atom in atoms:
        ans[atom.residue_id].append(atom)
    # append atoms {N, CA, C} of next residue
    for atom in atoms:
        if atom.name == "N" or atom.name == "CA" or atom.name == "C":
            chain, seq = atom.residue_id[0], int(atom.residue_id[1:])
            prev_residue_id = chain + str(seq - 1)
            if prev_residue_id in ans:
                newatom = atom._replace(residue_id = prev_residue_id, name = atom.name + "_next")
                ans[prev_residue_id].append(newatom)
    return ans


def get_residue_examples(atoms):
    residue_blks = get_residue_blocks(atoms)
    residue_blks = {
        residue_id: normalize_positions(residue_blks[residue_id])
        for residue_id in REPS_LIST
        if residue_id in residue_blks }
    ans = defaultdict(lambda: [])
    for residue_id in residue_blks:
        blk = residue_blks[residue_id]
        print(blk[0].residue_id, end="")
        ans[blk[0].residue].append(ResidueExample(
            [atom for atom in blk if atom.name not in {"CA_next", "C_next"}],
            get_transform(blk)))
    return ans


def atom_repr(atom):
    return "Atom(\"%s\", np.%s, \"%s\")" % (atom.name, repr(atom.pos), atom.element)

def residue_example_repr(eg):
    return "ResidueExample([\n            %s],\n        %s)" % (
        ",\n            ".join([atom_repr(atom) for atom in eg.atoms]),
        "np." + repr(eg.transform))

def main():
    atoms = parse_pdb(PDB_PATH)
    residue_examples = get_residue_examples(atoms)
    # clip so we don't end up with too many examples:
    for residue in residue_examples:
        residue_examples[residue] = residue_examples[residue][:10]
        print(residue, len(residue_examples[residue]))
    # code generation:
    print("\n\t--- Residue Construction Info: ---\n")
    print("{")
    for residue in residue_examples:
        print("    \"%s\": [" % residue)
        for eg in residue_examples[residue]:
            print("        %s," % residue_example_repr(eg))
        print("    ],")
    print("}")
    print()


if __name__ == "__main__":
    main()




