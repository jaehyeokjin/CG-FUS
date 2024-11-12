# Make PDB files for unfolded chains given an amino acid sequence

Start from an amino acid sequence and get a corresponding PDB file containing an unfolded protein with that sequence. The program makes an effort to avoid self-collisions and to generate a fairly compact chain. You can view the output and rerun if it is not satisfactory. Running an energy minimization afterwards is reccommended!

## Usage:

Pick a destination file. eg. `chains/8P59_unfolded.pdb`. Then run the program:

```
python seq2pdbchain.py > chains/8P59_unfolded.pdb
```

the program will wait for you to paste in an amino acid sequence (use the standard single-letter code) eg:

```
GSHMVPISFVFNRFPRMVRDLAKKMNKEVNFIMRGEDTELDRTFVEEIGEPLLHLLRNAIDHGIEPKEERIAKGKPPIGTLILSARHEGNNVVIEVEDDGRGIDKEKIIRKAIEKGLIDESKAATLSDQEILNFLFVPGFSTKEKVSEVSGRGVGMDVVKNVVESLNGSISIESEKDKGTKVTIRLPLT
```

Now it will run. Progress will be printed to stderr while the actual PDB data will be written to stdout. Be warned, this can take a while. Longer sequences tend to take more time than short ones.

```
reached 20 residues
reached 40 residues
reached 60 residues
reached 80 residues
reached 100 residues
reached 120 residues
reached 140 residues
reached 160 residues
reached 180 residues
generated radii: 72.66125837695655,  82.63373094440966,  74.6990087257651,  93.86458597438353,  82.11500262414255,  117.79293476432316,  93.02286714315417,  104.19893614030393,  60.73629381203191,  53.44501330294728,  102.62828107208239,  91.5770164819002,  92.72031625730851,  68.278014329689,  131.72060255164962,  61.492131423541316
best: 53.445013
```




