# rmsd.py
# Version 4.3
# Written by Kelsie M. King
# Github: kelsieking23
# Contact for issues and requests: kelsieking23@vt.edu
# Last updated: 01/06/2022
# Changes:
# * better help output
# * install pandas now optional
import argparse
from math import sqrt
import os
import sys
try:
    import pandas as pd
except:
    response = input('NOTE: The python package pandas is required for this script and is not installed. Download and install pandas? (type yes/no) ').lower()
    if (response.startswith('y')):
        print('\n')
        import pip
        pip.main(['install', 'pandas']) 
        import pandas as pd
    else:
        print('This script will now terminate. To use this script, install pandas using pip:\npip install pandas')

#TODO: maybe do clustering from RMSD matrix?
#TODO: write README & actual documentation, reference this in help

def getPoseCoords(ligand, split=False):
    '''
    Split poses in docked ligand .pdbqt file
    Arguments:
    * ligand (str): docking output ligand .pdbqt file
    * split (bool): whether to split ligand pdbqt into separate .pdbqts by model
    Returns: (dict {str:list}) dictionary mapping model name to a list of lists containing xyz coordiantes for model
    '''
    pose_lines = []
    pose_coords = []
    all_coords = {}
    all_lines = {}
    f = open(ligand, 'r')
    for line in f:
        line_parts = line.split()
        if ('HETATM' in line_parts[0]) or ('ATOM' in line_parts[0]):
            if line_parts[4].isalpha():
                try:
                    coords = list(map(float, line_parts[6:9]))
                except:
                    coords = fixBadCoordinates(line_parts[6:9])
            else:
                try:
                    coords = list(map(float, line_parts[5:8]))
                except:
                    coords = fixBadCoordinates(line_parts[5:8])
            pose_coords.append(coords)
            pose_lines.append(line)
        if 'MODEL' in line:
            current_model = line.strip()
        if 'ENDMDL' in line:
            all_coords[current_model] = pose_coords
            all_lines[current_model] = pose_lines
            pose_coords = []
            pose_lines = []
    if all_coords == {}:
        all_coords['MODEL 1'] = pose_coords
    f.close()
    return all_coords

def splitPoses(ligand):
    '''
    Write split poses into individual .pdbqt files.
    Arguments:
    * ligand (str): docked ligands .pdbqt file
    * coordinates (dict): coordinates by model
    Returns: (list) a list of output filenames
    '''
    filepaths = []
    current_model_lines = []
    f = open(ligand, 'r')
    for line in f:
        if 'MODEL' in line:
            current_model = ''.join(line.strip().split())
        elif 'ENDMDL' in line:
            filename = '{}_{}.pdbqt'.format(os.path.basename(ligand).split('.')[0], ''.join(current_model.split()))
            filepath = os.path.join(os.path.dirname(ligand), filename)
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
    f = open(reference, 'r')
    for line in f:
        line_parts = line.split()
        if ('HETATM' in line_parts[0]) or ('ATOM' in line_parts[0]):
            if line_parts[4].isalpha():
                try:
                    coords = list(map(float, line_parts[6:9]))
                except:
                    coords = fixBadCoordinates(line_parts[6:9])
            else:
                try:
                    coords = list(map(float, line_parts[5:8]))
                except:
                    coords = fixBadCoordinates(line_parts[5:8])
            ref_coords.append(coords)
    f.close()
    return ref_coords

def fixBadCoordinates(line_parts):
    '''
    Fixes bad coordinates with no spaces between them
    Arguments:
    * line_parts (list): list of bad coordinates
    Returns: (list) fixed x,y,z coordinates for atom
    '''
    coords = []
    chunk = ''.join(line_parts)
    string = ''
    i = 0
    for char in chunk:
        if (char != '.') and (i == 0):
            string = string + char
        if (char != '.') and (i != 0):
            string = string + char
            i += 1
            if i == 4:
                coords.append(float(string))
                if len(coords) == 3:
                    return coords
                i = 0
                string = ''
        if char == '.':
            string = string + char
            i += 1

def rmsd(ligand, reference, output=None):
    '''
    Calculates RMSD
    Arguments:
    * ligand (str): ligand .pdbqt file
    * reference (str): reference .pdbqt file
    * output (optional, str): output filepath (.dat or .csv)
    Returns: (list) list of tuples (model name, rmsd)
    '''
    rmsds = []
    all_ligand_coords = getPoseCoords(ligand)
    ref_coords = getRefCoords(reference)
    for key in all_ligand_coords.keys():
        ligand_coords = all_ligand_coords[key]
        len_coords = len(ligand_coords)
        total_squared_distance = 0
        for i in range(0, len_coords):
            total_squared_distance += (ligand_coords[i][0] - ref_coords[i][0])**2 + (ligand_coords[i][1] - ref_coords[i][1])**2 + (ligand_coords[i][2] - ref_coords[i][2])**2
        rmsd = sqrt(total_squared_distance / len_coords)
        data = (key, rmsd)
        rmsds.append(data)
    if output is not None:
        ind = []
        vals = []
        for item in rmsds:
            ind.append(''.join(item[0].split()))
            vals.append(item[1])
        df = pd.DataFrame()
        df['rmsd'] = vals
        df.index=ind
        df.to_csv(output)
    return rmsds

def rmsdMatrix(ligand, output=None):
    '''
    Outputs RMSD matrix between all ligand poses.
    Arguments:
    * ligand (str): ligand .pdbqt file
    * output (optional, str): output filepath (.csv)
    Returns: (pandas DataFrame) RMSD matrix
    '''
    rmsds = {}
    all_ligand_coords = getPoseCoords(ligand)
    for key1 in all_ligand_coords.keys():
        ligand_coords = all_ligand_coords[key1]
        if key1 not in rmsds.keys():
            rmsds[key1] = {}
        for key2 in all_ligand_coords.keys():
            ref_coords = all_ligand_coords[key2]
            len_coords = len(ligand_coords)
            total_squared_distance = 0
            for i in range(0, len_coords):
                total_squared_distance += (ligand_coords[i][0] - ref_coords[i][0])**2 + (ligand_coords[i][1] - ref_coords[i][1])**2 + (ligand_coords[i][2] - ref_coords[i][2])**2
            rmsd = sqrt(total_squared_distance / len_coords)
            rmsds[key1][key2] = rmsd
    df = pd.DataFrame(rmsds)
    if output is not None:
        df.to_csv(output)
    return df
    
def errorHandler(ligand, reference, args):
    '''
    Error handling
    '''
    # check if .pdbqt
    if ligand is not None:
        if not ligand.endswith('pdbqt'):
            raise ValueError("Ligand file {} is not a .pdbqt. Please check inputs.".format(ligand))
        with open(ligand, 'r') as f:
            ligand_contents = f.readlines()
        if ligand_contents == []:
            raise ValueError('Ligand file {} is empty. Please check inputs.')
    else:
        raise ValueError('Ligand (-l) must be specified.')

    if reference is not None:
        if not reference.endswith('pdbqt'):
            raise ValueError("Reference file {} is not a .pdbqt. Please check inputs.".format(reference))
        with open(reference, 'r') as f:
            reference_contents = f.readlines()
        if reference_contents == []:
            raise ValueError('Reference file {} is empty. Please check inputs.')

    # check if output filetype is valid
    if (args.__dict__['matrix'] is not 'None') and (args.__dict__['matrix'] is not None):
        if (not args.__dict__['matrix'].endswith('.csv')) and (not args.__dict__['matrix'].endswith('.dat')):
            raise ValueError('Output file for matrix ({}) is not a valid file type (must be .csv or .dat)'.format(args.__dict__['matrix']))
    if args.__dict__['o'] is not None:
        if (not args.__dict__['o'].endswith('.csv')) and (not args.__dict__['o'].endswith('.dat')):
            raise ValueError('Output file for matrix ({}) is not a valid file type (must be .csv or .dat)'.format(args.__dict__['o']))

    # check that number of atoms are the same
    if (ligand is not None) and (reference is not None):
        ligand_atoms = 0
        reference_atoms = 0
        for line in ligand_contents:
            line_parts = line.split()
            if ('HETATM' in line_parts[0]) or ('ATOM' in line_parts[0]):
                ligand_atoms += 1
            if 'ENDMDL' in line:
                break
        for line in reference_contents:
            line_parts = line.split()
            if ('HETATM' in line_parts[0]) or ('ATOM' in line_parts[0]):
                reference_atoms += 1
        if ligand_atoms != reference_atoms:
            raise ValueError('Number of atoms in ligand file ({}) do not match number of atoms in reference file ({}). Please check inputs.'.format(ligand_atoms, reference_atoms))

if __name__ == '__main__':
    # create parser
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
    parser.add_argument('-l', help='(REQUIRED) Path to docking output file (.pdbqt). This is the file containing all docked poses.')
    parser.add_argument('-r', help='Path to reference pose file (.pdbqt) This is the original ligand inputted to Vina.')
    parser.add_argument('-o', nargs='?', help='Output file (.dat or .csv)')
    parser.add_argument('--matrix', '-m', nargs='?', default='None',
                        help='If specified, output a matrix of RMSD values between all poses in the docking output file. Can specify a .csv or .dat filename to save matrix to a file')
    parser.add_argument('--split', '-s', action='store_true', help='If specified, will split poses in docking output (-l) into individual .pdbqt files')
    args = parser.parse_args()


    # show help message if ran with no arguments
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # parse args & error handle input files
    ligand = args.__dict__['l']
    reference = args.__dict__['r']
    errorHandler(ligand, reference, args)

    # split files if specified
    if args.__dict__['split'] is True:
        filepaths = splitPoses(ligand)
        print('Wrote {} files:'.format(len(filepaths)))
        for filepath in filepaths:
            print('  {}'.format(filepath))

    # do rmsd
    if (ligand is not None) and (reference is not None):
        if args.__dict__['o'] is not None:  
            output = args.__dict__['o']
            rmsds = rmsd(ligand, reference, output=output)
        else:
            output = None
            rmsds = rmsd(ligand, reference)
        print('********')
        print('RMSD compared to reference ligand {}'.format(os.path.abspath(reference)))
        if output is not None:
            print('Saved RMSD to {}'.format(os.path.abspath(output)))
        print('********')
        for item in rmsds:
            print(item[0] + ' RMSD: ' + str(item[1]))

    # do RMSD matrix
    if args.__dict__['matrix'] is not 'None':
        if args.__dict__['matrix'] is None:
            matrix_output = os.path.join(os.getcwd(), 'rmsd_matrix.csv')
            df = rmsdMatrix(ligand, output=matrix_output)
        else:
            matrix_output = args.__dict__['matrix']
            df = rmsdMatrix(ligand, output=matrix_output)
        print('********')
        print('RMSD matrix')
        print('Saved to {}'.format(os.path.abspath(matrix_output)))
        print('********')
        print(df)