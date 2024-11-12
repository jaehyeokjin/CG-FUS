# Step 1: Atomistic Setup for FUS Configuration

# 1. Codes
- The code files in `./codes` generate the all-atomistic (united-atom level; hydrogen atoms are omitted) and coarse-grained (carbon alpha-level) structure of a general protein.
- The base code is adapted from [seq2pdbchain](https://github.com/pb1729/seq2pdbchain) to avoid self-collisions between amino acids while generating a compact amino acid chain.
- This code produces an unfolded chain, making it suitable for studying intrinsically disordered proteins (IDP) and intrinsically disordered regions (IDR). Note that energy minimization is still recommended (which we will perform in Step 2 in the subsequent folder) for single-chain level descriptions.
- Several fixes and typos from the original GitHub repository were corrected. The forked version can be found [here](link-to-forked-repo).

## 2. FUS Protein (How to Use)
- To generate AA and CG configurations of the FUS protein, run the following command:
  `python seq2pdbchain.py`
- When prompted, provide the following amino acid sequence:
	MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNS
- The script will generate the following PDB files:
  - `./fus.pdb` for the all-atom (AA) structure.
  - `./cg_fus.pdb` for the coarse-grained (CG) structure.
- Use PDB rendering software to visualize the generated PDB files.

## 3. Miscellaneous
- If needed, `convert.py` allows for converting `fus.pdb` (AA level) to `*.xyz` format for further analysis. An example output, `output.xyz`, is provided.

