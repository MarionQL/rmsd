# rmsd
Performs RMSD calculations for redocking.

**rmsd.py** takes the outputted structure file from docking with Vina (**-l**) and the original ligand file, or reference file (**-r**) and outputs RMSD values. The outputted RMSD values can be saved to a file by specifying an output path (**-o**). This script can also output an RMSD matrix between all poses in the outputted structue file (**-m**). Additionally, the outputted structure file from vina (**-l**) can be split into indivindual .pdbqt files using the option **-s**.

## Usage
**-l** _(.pdbqt)_ Path to docking output file. This is the file containing all docked poses.

**-r** _(.pdbqt)_ Path to reference structure file. This is the original structure file inputted to Vina.

**-o** _(.csv, .dat, default = None)_ Output file for RMSD values.

**-s, --split** _(default = False)_ If specified, will split poses in docking output file (specified with **-l**) into individual .pdbqt files.

**-m, --matrix** _(.csv, .dat, default = None)_ If specified, output a matrix of RMSD values between all poses in the docking output file (specified with **-l**). Can specify an output file to save matrix. 

**-h, --help** Show help message and quit.


## Examples
To compute the RMSD between output from Vina and the reference structure:
```
>>> python rmsd.py -l docking_output.pdbqt -r reference_file.pdbqt -o rmsd.csv
********
RMSD compared to reference ligand docking_output.pdbqt
Saved RMSD to rmsd.csv
********
MODEL 1 RMSD: 2.704676149573933
MODEL 2 RMSD: 2.7385713699204826
MODEL 3 RMSD: 6.420578703581858
MODEL 4 RMSD: 4.49323453567612
MODEL 5 RMSD: 3.76852366711815
MODEL 6 RMSD: 3.948213248021421
MODEL 7 RMSD: 3.2767659304507792
MODEL 8 RMSD: 4.031011333446035
MODEL 9 RMSD: 6.698428829044747
```

To compute the RMSD between output from Vina and the reference structure, and also split the docked poses file:
```
>>> python rmsd.py -l docking_output.pdbqt -r reference_file.pdbqt -o rmsd.csv --split
Wrote 9 files:
  docking_output_MODEL1.pdbqt
  docking_output_MODEL2.pdbqt
  docking_output_MODEL3.pdbqt
  docking_output_MODEL4.pdbqt
  docking_output_MODEL5.pdbqt
  docking_output_MODEL6.pdbqt
  docking_output_MODEL7.pdbqt
  docking_output_MODEL8.pdbqt
  docking_output_MODEL9.pdbqt
********
RMSD compared to reference ligand reference_file.pdbqt
Saved RMSD to rmsd.csv
********
MODEL 1 RMSD: 2.704676149573933
MODEL 2 RMSD: 2.7385713699204826
MODEL 3 RMSD: 6.420578703581858
MODEL 4 RMSD: 4.49323453567612
MODEL 5 RMSD: 3.76852366711815
MODEL 6 RMSD: 3.948213248021421
MODEL 7 RMSD: 3.2767659304507792
MODEL 8 RMSD: 4.031011333446035
MODEL 9 RMSD: 6.698428829044747
```
