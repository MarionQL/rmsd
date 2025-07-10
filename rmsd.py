# rmsd.py
# Version 4.7
# Written by Kelsie M. King
# Github: kelsieking23
# Contact for issues and requests: kelsieking23@vt.edu
# Last updated: 7/10/2025
# Changes:
# * bug fix for output file
# Symmetric RMSD calculations and refactoring by Marion LoPresti

import argparse
from math import sqrt
import os
import sys
import subprocess
import warnings

warnings.filterwarnings("ignore")

# Auto-installation helpers
def install_package(package_name, install_cmd, submodules=None, manager_name=None):
    '''
    Attempts to import a Python package, and installs it if missing.

    Arguments:
    * package_name (str): The name of the package to import.
    * install_cmd (list): The command to install the package (e.g., ['pip', 'install', 'package']).
    * submodules (list of str, optional): List of submodules to import after installation,
      optionally with aliases (e.g., ['pandas as pd']).
    * manager_name (str, optional): Name of the package manager for user prompts (e.g., 'pip' or 'conda').

    Returns:
    * modules (dict): A dictionary mapping submodule names (or aliases) to imported modules.

    Notes:
    - Prompts the user before installing missing packages.
    - Displays specific guidance for RDKit or Brown Lab students.
    - Automatically imports specified submodules after installation.
    '''
    try:
        __import__(package_name)
    except ImportError:
        response = input(f'{package_name} not found. Install it via {manager_name or "pip"}? (yes/no) ').lower()
        if response.startswith('y'):
            try:
                subprocess.check_call(install_cmd)
            except Exception:
                print(f'Error: Installation of {package_name} failed.\n'
                      f'Please try installing it manually with the following command:\n\n'
                      f'    {" ".join(install_cmd)}\n')
                
                # Brown Lab–specific note for RDKit or Conda
                if "conda" in install_cmd[0] or package_name.lower() == 'rdkit':
                    print(
                        '\n***If you are a Brown Lab student and manual installation fails:\n'
                        '\tThis script requires an Anaconda installation and Python 3 on your computer to work.\n'
                        '\t> On the Brown Lab canvas page, visit BL Module 3 | Indroduction to Computational Literacy\n'
                        '\t> On page (3/4), follow the instructions for installing Anaconda Navigator on your computer\n'
                        '\t> After installation, close your current terminal session and try running this script again.\n'
                        '***If you are NOT a Brown Lab student and manual installation fails:\n'
                        '\t> Download and install Anaconda (https://www.anaconda.com/products/distribution)\n'
                        '\t> Close your current terminal session and try running this script again.'
                    )
                sys.exit(1)
        else:
            print(f'Please install {package_name} manually to continue.')
            sys.exit(1)
    
    modules = {}
    if submodules:
        for submodule in submodules:
            parts = submodule.split(" as ")
            import_stmt = parts[0].strip()
            alias = parts[1].strip() if len(parts) == 2 else None
            
            module = __import__(import_stmt, fromlist=[''])
            modules[alias or import_stmt] = module
    return modules

modules = install_package(
    'pandas',
    [sys.executable, "-m", "pip", "install", "pandas"],
    submodules=['pandas as pd'],
    manager_name='pip'
)

pd = modules['pd']  # Just like `import pandas as pd`

# RDKit
rdkit_modules = install_package(
    'rdkit',
    ['conda', 'install', '-c', 'conda-forge', 'rdkit', '-y'],
    submodules=[
        'rdkit.Chem as Chem',
        'rdkit.Chem.rdMolAlign as rdMolAlign',
        'rdkit.Geometry as Geometry'
    ],
    manager_name='conda'
)

Chem = rdkit_modules['Chem']
rdMolAlign = rdkit_modules['rdMolAlign']
Point3D = rdkit_modules['Geometry'].Point3D

def parse_coords(line):
    '''
    Extracts the x, y, z coordinates from a PDB-format line.
    Arguments:
    * line (str): A line from a PDB or PDBQT file starting with 'ATOM' or 'HETATM'.
    Returns:
    * coords (list of float): A list containing the x, y, z coordinates as floats.
    Notes:
    - Assumes standard PDB column positions for coordinates (columns 31–54).
    '''
    try:
        return [float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())]
    except Exception as e:
        raise ValueError(f"Failed to parse coordinates: {line}\nError: {e}")

def getPoseCoords(ligand):
    '''
    Parses a ligand file containing multiple poses (e.g., docking output in .pdbqt or .pdb format)
    and extracts the 3D coordinates for each pose.
    Arguments:
    * ligand (str): Path to the ligand file (.pdbqt or .pdb) containing multiple poses.
                    The file should contain standard MODEL/ENDMDL blocks or equivalent headers
                    to delimit poses (e.g., MODEL, HEADER, CRYST1, ENDMDL, END).
    Returns:
    * poses (dict): Dictionary mapping model names (e.g., 'MODEL_1', 'MODEL_2', ...) to lists of
                    atomic coordinates. Each value is a list of [x, y, z] coordinates for all atoms
                    in that pose.
    Notes:
    - Atom coordinates are parsed from 'HETATM' and 'ATOM' lines.
    - Poses are auto-numbered if MODEL records are missing.
    - Handles both multi-pose and single-pose files.
    '''
    poses = {}
    coords = []
    pose_idx = 1
    current_model = None

    with open(ligand, 'r') as f:
        for line in f:
            if any(tag in line for tag in ['MODEL', 'HEADER', 'CRYST1']):
                if coords:
                    poses[current_model or f'MODEL_{pose_idx}'] = coords
                    coords = []
                    pose_idx += 1
                current_model = f'MODEL_{pose_idx}'
            elif line.startswith(('HETATM', 'ATOM')):
                coords.append(parse_coords(line))
            elif 'ENDMDL' in line or 'END' in line:
                if coords:
                    poses[current_model or f'MODEL_{pose_idx}'] = coords
                    coords = []
                    pose_idx += 1
                    current_model = None

    if coords:
        poses[current_model or f'MODEL_{pose_idx}'] = coords
    return poses

def splitPoses(ligand):
    '''
    Write split poses into individual .pdbqt files.
    Arguments:
    * ligand (str): docked ligands .pdbqt file
    Returns: (list) a list of output filenames
    '''
    filepaths = []
    current_model_lines = []
    f = open(ligand, 'r')
    for line in f:
        if 'MODEL' in line or 'HEADER' in line or 'CRYST1' in line:
            current_model = ''.join(line.strip().split())
        elif 'ENDMDL' in line or 'END' in line:
            base, ext = os.path.splitext(ligand)
            filename = '{}_{}{}'.format(base, ''.join(current_model.split()), ext)
            output_dir = os.path.dirname(ligand) or '.'
            filepath = os.path.join(output_dir, filename)
            filepaths.append(filepath)
            with open(filepath, 'w') as g:
                for line in current_model_lines:
                    g.write(line)
            current_model_lines = []
        else:
            current_model_lines.append(line)
    f.close()
    return filepaths

def getRefCoords(reference):
    '''
    Gets reference coords from .pdbqt file
    Arguments: 
    * reference (str): reference ligand .pdbqt file
    Returns: (list) list of lists containing xyz coordinates for reference 
    '''
    ref_coords = []
    with open(reference, 'r') as f:
        for line in f:
            if line.startswith(('HETATM', 'ATOM')):
                ref_coords.append(parse_coords(line))
    return ref_coords

def rmsd(ligand, reference):
    '''
    Calculates RMSD
    Arguments:
    * ligand (str): ligand .pdbqt file
    * reference (str): reference .pdbqt file
    Returns: (list) list of tuples (model name, rmsd)
    '''
    ligand_coords = getPoseCoords(ligand)
    ref_coords = getRefCoords(reference)
    rmsds = []

    for model, coords in ligand_coords.items():
        total = sum((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2 for (x1, y1, z1), (x2, y2, z2) in zip(coords, ref_coords))
        rmsds.append((model, sqrt(total / len(coords))))
    return rmsds

def build_rdkit_mol(coords):
    '''
    Builds an RDKit molecule object with a single conformer from a list of 3D coordinates.
    Arguments:
    * coords (list of lists): List of [x, y, z] coordinate triplets for each atom.
    Returns:
    * RDKit Mol object with dummy carbon atoms (atomic number 6) at the specified coordinates.
    Notes:
    - All atoms are treated as carbons; only 3D coordinates are encoded.
    - Useful for RMSD calculations where atom types are irrelevant.
    '''
    mol = Chem.RWMol()
    conf = Chem.Conformer(len(coords))
    for idx, (x, y, z) in enumerate(coords):
        mol.AddAtom(Chem.Atom(6))
        conf.SetAtomPosition(idx, Point3D(x, y, z))
    mol.AddConformer(conf)
    return mol.GetMol()

def rdkit_rmsd(ligand, reference):
    '''
    Calculates RMSD between each ligand pose and a reference structure using RDKit's CalcRMS,
    which accounts for molecular symmetry but does not modify coordinates.
    Arguments:
    * ligand (str): Path to ligand file (.pdbqt or .pdb) containing multiple poses.
    * reference (str): Path to reference ligand file (.pdbqt or .pdb).
    Returns:
    * list of tuples: [(pose_name, RMSD_value)] for each pose in the ligand file.
    Notes:
    - Uses RDKit’s CalcRMS, which considers all symmetry-equivalent mappings.
    - Molecule coordinates are not altered during calculation.
    '''
    ligand_coords = getPoseCoords(ligand)
    ref_mol = build_rdkit_mol(getRefCoords(reference))
    return [(model, rdMolAlign.CalcRMS(build_rdkit_mol(coords), ref_mol)) for model, coords in ligand_coords.items()]

def writeRMSDToFile(rmsds, output):
    '''
    Writes RMSD values to a file. File is written in .csv format (model, rmsd)
    Arguments:
    * rmsds (list): Output of rmsd(): list of tuples [(model name, rmsd)]
    * output (str): output filepath 
    '''
    with open(output, 'w') as f:
        f.write('model,rmsd\n')
        for model, rmsd_val in rmsds:
            f.write(f"{model},{rmsd_val:.3f}\n")

def rmsdMatrix(ligand, use_rdkit=False, output=None):
    '''
    Outputs RMSD matrix between all ligand poses.
    Arguments:
    * ligand (str): ligand .pdbqt file
    * use_rdkit (bool): If True, the RDKit symmetry aware RMSD will be used
    * output (optional, str): output filepath (.csv)
    Returns: (pandas DataFrame) RMSD matrix
    '''
    ligand_coords = getPoseCoords(ligand)
    models = list(ligand_coords.keys())
    matrix = pd.DataFrame(index=models, columns=models, dtype=float)

    for i, model1 in enumerate(models):
        for j, model2 in enumerate(models):
            if use_rdkit:
                mol1 = build_rdkit_mol(ligand_coords[model1])
                mol2 = build_rdkit_mol(ligand_coords[model2])
                rmsd_val = rdMolAlign.CalcRMS(mol1, mol2)
            else:
                total = sum((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2 for (x1, y1, z1), (x2, y2, z2) in zip(ligand_coords[model1], ligand_coords[model2]))
                rmsd_val = sqrt(total / len(ligand_coords[model1]))
            
            matrix.at[model1, model2] = rmsd_val

    if output:
        matrix.to_csv(output)
    return matrix

def errorHandler(ligand, reference, args):
    '''
    Error handling
    '''
    # check if .pdbqt or .pdb
    if ligand is not None:
        if not ligand.endswith(('pdbqt', 'pdb')):
            raise ValueError(f"Ligand file {ligand} is not a .pdbqt or .pdb. Please check inputs.")
        with open(ligand, 'r') as f:
            ligand_contents = f.readlines()
        if not ligand_contents:
            raise ValueError(f'Ligand file {ligand} is empty. Please check inputs.')
    else:
        raise ValueError('Ligand (-l) must be specified.')

    if reference is not None:
        if not reference.endswith(('pdbqt', 'pdb')):
            raise ValueError(f"Reference file {reference} is not a .pdbqt or .pdb. Please check inputs.")
        with open(reference, 'r') as f:
            reference_contents = f.readlines()
        if not reference_contents:
            raise ValueError(f'Reference file {reference} is empty. Please check inputs.')

    # check that number of atoms are the same
    if ligand and reference:
        all_ligand_coords = getPoseCoords(ligand)
        ligand_atom_counts = {pose: len(coords) for pose, coords in all_ligand_coords.items()}
        ref_coords = getRefCoords(reference)
        reference_atoms = len(ref_coords)
        
        for pose_name, atom_count in ligand_atom_counts.items():
            if atom_count != reference_atoms:
                raise ValueError(
                    f'Number of atoms in ligand pose "{pose_name}" ({atom_count}) '
                    f'does not match reference ({reference_atoms}). Please check inputs.'
                )

    # Check output directory exists
    if args.o:
        directory = os.path.dirname(args.o)
        if directory and not os.path.isdir(directory):
            raise ValueError(f'No such directory: {directory}\nPlease check output path.')

def main():
    description = '''Performs RMSD calculations.
        * Compute RMSD between reference ligand and docking outputs for redocking
        * Compute RMSD matrix between all docked poses in a .pdbqt file
Inputs must be .pdbqt'''
    epilog =  '''Examples:
Example redocking command: python rmsd.py -l docking_output.pdbqt -r reference_pose.pdbqt
Example redocking command, saving output to file: python rmsd.py -l docking_output.pdbqt -r reference_pose.pdbqt -o output.csv
Example RMSD matrix command: python rmsd.py -l docking_output.pdbqt -m rmsd_matrix.csv
Example RMSD matrix + redocking command: python rmsd.py -l docking_output.pdbqt -r reference_pose.pdbqt -m rmsd_matrix.csv
Example RMSD matrix + redocking command, saving redocking output: python rmsd.py -l docking_output.pdbqt -r reference_pose.pdbqt -m rmsd_matrix.csv -o output.csv'''
    parser = argparse.ArgumentParser(description=description, epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-l', help='(REQUIRED) Path to docking output file (.pdbqt or .pdb). This is the file containing all docked poses.')
    parser.add_argument('-r', help='Path to reference pose file (.pdbqt or .pdb) This is the original ligand inputted to Vina.')
    parser.add_argument('-o', nargs='?', help='Path to output file for RMSD values')
    parser.add_argument('--matrix', '-m', nargs='?', const=True, default=None,
                        help='If specified, output a matrix of RMSD values between all poses in the docking output file. Can specify a .csv or .dat filename to save matrix to a file')
    parser.add_argument('--split', '-s', action='store_true', help='If specified, will split poses in docking output (-l) into individual .pdbqt files')
    parser.add_argument('--rdkit', '-k', action='store_true', help='Use RDKit-based permutation-invariant RMSD instead of simple coordinate-based RMSD.')

    args = parser.parse_args()

    errorHandler(args.l, args.r, args)
    if args.r:
        print('******** RMSD Calculation ********')
        method = 'RDKit-based RMSD (symmetry-aware)' if args.rdkit else 'Simple coordinate RMSD'
        print(f'RMSD Method: {method}')
        if args.rdkit:
            rmsds = rdkit_rmsd(args.l, args.r)
        else:
            rmsds = rmsd(args.l, args.r)

        if args.o:
            writeRMSDToFile(rmsds, args.o)
            print(f'\nRMSD values saved to: {os.path.abspath(args.o)}')
        else:
            for model, val in rmsds:
                print(f"{model} RMSD: {val:.3f}")
    
    if args.split:
        print('******** Pose Splitting ********')
        for filepath in filepaths:
            print(f'  - {filepath}')
        filepaths = splitPoses(args.l)
        print('Wrote {} files:'.format(len(filepaths)))
        for filepath in filepaths:
            print('  {}'.format(filepath))

    if args.matrix:
        output_file = None if args.matrix is True else args.matrix
        matrix = rmsdMatrix(args.l, use_rdkit=args.rdkit, output=output_file)
        print("\n******** RMSD Matrix Calculation ********")
        print("\nRMSD Matrix:")
        print(matrix)
        if isinstance(args.matrix, str):
            print(f'\nMatrix saved to: {os.path.abspath(args.matrix)}')
        else:
            print('\nMatrix not saved (no filename specified).')

if __name__ == '__main__':
    main()
