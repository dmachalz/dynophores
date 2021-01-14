"""This is a helper script to generate renaming dictionaries for dynophores analysis"""

import os
import sys
import configparser
import pickle

sys.path.extend(['/mdspace/davidm/tools/dynophores/analysis'])  # has to be changed, if path changes
from lib.functions import getligatomdict
from lib.functions import getresatomdict

# Reading in config file
configpath = '/home/davidm/Tools/dynophores/testing/test-getligrenamedicts.cfg'
if len(sys.argv) == 2 and os.path.exists(sys.argv[1]):
    configpath = sys.argv[1]  # Config file path provided via command line
    print('Using config file {}'.format(sys.argv[1]))
elif configpath is not None:
    pass  # this is used for the case, that the configfile is defined within the script
else:
    config = None
    print('No (suitable) config file provided!')
    exit(1)

config = configparser.ConfigParser()
config.sections()
config.read(configpath)

# Interpreting config
if config is not None and config.has_section('Renaming dictionaries'):

    if config.get('Renaming dictionaries', 'list of pdb data paths') != '':
        # Input should be tested - like raw data paths
        pdbinlist = config.get('Renaming dictionaries', 'list of pdb data paths').split(' ')
    else:
        pdbinlist = None
        print('No pdb files found')
        exit()

    if config.get('Renaming dictionaries', 'datalabels') != '':
        datalabels = config.get('Renaming dictionaries', 'datalabels').split(' ')
    else:
        datalabels = None  # This will tell the readdynos function to generate generic datalabels raw_1...
        print('Warning: no datalabels defined')
        exit()

    if config.has_option('Renaming dictionaries', 'ligand name tag') and \
            config.has_option('Renaming dictionaries', 'ligand rename dictionary path'):
        if config.get('Renaming dictionaries', 'ligand name tag') != '' and \
                config.get('Renaming dictionaries', 'ligand rename dictionary path') != '':
            ligname = config.get('Renaming dictionaries', 'ligand name tag')
            ligatomdictpath = config.get('Renaming dictionaries', 'ligand rename dictionary path')
        else:
            ligname = ligatomdictpath = None
    else:
        ligname = ligatomdictpath =None

    if config.has_option('Renaming dictionaries', 'residues rename dictionary path'):
        resatomdictpath = config.get('Renaming dictionaries', 'residues rename dictionary path')
    else:
        resatomdictpath = None

else:
    pdbinlist = datalabels = ligname = ligatomdictpath = resatomdictpath = None
    print('No suitable config file found')
    exit()

# Reading in files
pdbtuplelist = []
if len(pdbinlist) == len(datalabels):
    for x, y in zip(pdbinlist, datalabels):
        pdbtuplelist.append((x, y))
else:
    print('Length of pdb file list is not the same length as datalabel list')
    exit()

# Getting dictionaries
if ligatomdictpath is not None and ligname is not None:
    print('Getting ligand atom dictionary')
    ligdict = getligatomdict(pdbtuplelist, ligname)
    # Checking for existing dictionary
    if os.path.exists(ligatomdictpath):
        print('Appending to ligand atom dictionary under {}'.format(ligatomdictpath))
        oldligdict = pickle.load(open(ligatomdictpath, "rb"))
        newligdict = {**oldligdict, **ligdict}  # Updating old dictionary with new entries
        pickle.dump(newligdict, open(ligatomdictpath, "wb"))  # Writing out new updated dictionary

    else:
        print('Writing out ligand atom dictionary to {}'.format(ligatomdictpath))
        pickle.dump(ligdict, open(ligatomdictpath, "wb"))

if resatomdictpath is not None:
    print('Getting residue atom dictionary')
    resdict = getresatomdict(pdbtuplelist)
    # Checking for existing dictionary
    if os.path.exists(resatomdictpath):
        print('Appending to residue atom dictionary under {}'.format(resatomdictpath))
        oldresdict = pickle.load(open(resatomdictpath, "rb"))
        newresdict = {**oldresdict, **resdict}  # Updating old dictionary with new entries
        pickle.dump(newresdict, open(resatomdictpath, "wb"))  # Writing out new updated dictionary
    else:
        print('Writing out ligand atom dictionary to {}'.format(resatomdictpath))
        pickle.dump(resdict, open(resatomdictpath, "wb"))
