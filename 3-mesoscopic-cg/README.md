# Step 3: Mesoscopic Coarse-Graining Mapping of FUS

### 1. General Explanation
- Building on molecular coarse-graining, mesoscopic coarse-graining further simplifies the coarse-grained particles into discretized fields.
- This approach involves applying Gaussian smearing to the microscopic Dirac delta operator. Specifically, we replace the microscopic Dirac delta operator localized on each particle, i.e., \(\rho(r) = \sum_i \delta(r - r_i)\), with a coarse-grained density function \(\rho_{CG}(r) := \sum_i G(r_i, \sigma)\), where \(G\) defines a Gaussian kernel centered at \(r_i\) with variance \(\sigma^2\). The parameter \(\sigma\) determines the level of mesoscopic coarse-graining.

### 2. Prerequisites and Codes
- This step requires the coarse-grained (CG) molecular dynamics (MD) simulation trajectory of FUS condensation from Step 2.
- We use the sampled trajectory from Step 2.
- Note that the trajectory in Step 2 was in `.dcd` format. For the purposes of analyzing the CG trajectory and performing mesoscopic coarse-graining, the LAMMPS trajectory format (`.lammpstrj`) is more convenient. Therefore, we converted `fus_120_run.dcd` to `final.lammpstrj`.
- These files are available in the Zenodo repository.

### 3. Usage
- Due to the large file size, this GitHub folder only includes a specific snapshot at the 933rd timestep. This file can be found as `933.lammpstrj`.
- To perform the mesoscopic coarse-graining, open and run the Jupyter Notebook file `Mesoscopic Coarse-Graining.ipynb`.
- This notebook takes the initial `933.lammpstrj` file, parses the coordinates, and then performs mesoscopic coarse-graining.
- To effectively smear out microscopic correlations (note that as \(\sigma \to 0\), it recovers the microscopic configuration), we chose \(\sigma_{md} = 15\).
- Additionally, we discretize the system into grid boxes. Here, we use a 50x500 grid (with a 1:10 aspect ratio based on the x/z length ratio) to discretize the system into a structured grid.


