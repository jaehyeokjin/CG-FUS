import numpy as np
from collections import namedtuple


PDBAtom = namedtuple("PDBAtom", ["name", "residue", "residue_id", "pos", "element"])
ResidueExample = namedtuple("ResidueExample", ["atoms", "transform"])


def get_record(line):
    """ https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html """
    if len(line) >= 78 and line[:4] == "ATOM":
        name = line[12:16].strip()
        residue = line[17:20].strip()
        residue_id = "%s%d" % (line[21], int(line[22:26]))
        chain_id = line[21]
        pos = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        element = line[76:78].strip()
        return PDBAtom(name, residue, residue_id, pos, element)
    return None


def parse_pdb(fnm):
    ans = []
    with open(fnm, "r") as f:
        for line in f.readlines():
            record = get_record(line)
            if record is not None:
                ans.append(record)
    return ans


def write_xyz(fnm, atoms):
    with open(fnm, "w") as f:
        f.write("%d\n\n" % len(atoms))
        for atom in atoms:
            f.write("%s   %f   %f   %f\n" % (atom.element, atom.pos[0], atom.pos[1], atom.pos[2]))


def leftpad(s, n):
    return " "*(n-len(s)) + s

def pdb_name(s):
    return " " + s + " "*(3-len(s))

def pdb_float(f):
    return "%8.3f" % f

def pdb_line(atom, i_atom, i_seq, residue):
    """ based on this: https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html """
    return ("ATOM  "                # 1-6
        + leftpad(str(i_atom), 5)   # 7-11
        + " "                       # 12
        + pdb_name(atom.name)       # 13-16
        + " "                       # 17
        + residue                   # 18-20
        + " A"                      # 21-22
        + leftpad(str(i_seq), 4)    # 23-26
        + "    "                    # 27-30
        + pdb_float(atom.pos[0])    # 31-38
        + pdb_float(atom.pos[1])    # 39-46
        + pdb_float(atom.pos[2])    # 47-54
        + "  1.00"                  # 55-60
        + "100.00"                  # 61-66
        + "          "              # 67-76
        + leftpad(atom.element, 2)  # 77-78
        + "  "                      # 79-80
    )


