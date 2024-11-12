def pdb_to_xyz(pdb_filename, xyz_filename):
    atoms = []
    
    with open(pdb_filename, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                element = line[76:78].strip()  # Element symbol
                x = float(line[30:38].strip())  # X coordinate
                y = float(line[38:46].strip())  # Y coordinate
                z = float(line[46:54].strip())  # Z coordinate
                atoms.append((element, x, y, z))
    
    with open(xyz_filename, 'w') as xyz_file:
        # Write the number of atoms
        xyz_file.write(f"{len(atoms)}\n")
        xyz_file.write(f"Converted from {pdb_filename}\n")
        # Write each atom's element and coordinates
        for atom in atoms:
            element, x, y, z = atom
            xyz_file.write(f"{element} {x:.3f} {y:.3f} {z:.3f}\n")

# Usage:
pdb_filename = 'fus.pdb'  # Replace with your PDB file
xyz_filename = 'output.xyz'  # Output XYZ file

pdb_to_xyz(pdb_filename, xyz_filename)

