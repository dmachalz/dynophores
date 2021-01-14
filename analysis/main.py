"""Analysing dynophore raw data in the desired way as specified in the config file"""

import os
import pickle
import sys
import ast
import configparser

sys.path.extend(['/mdspace/davidm/tools/dynophores/analysis'])  # has to be changed, if path changes
from lib.functions import *


# Explicit testing variables
figuresize = (5, 8)
configpath = '/mdspace/davidm/cyp4z1/mdsims/analysis/dyno-conffiles/dyno-analysis-wt-1-distances.cfg'

### Stage (0) ###
# Reading in config file
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

### Processing config file ###

if config is not None and config.has_section('Reading and Processing'):
    # (A) Raw data paths
    if config.has_option('Reading and Processing', 'List of raw data paths') and \
            config.get('Reading and Processing', 'List of raw data paths') != '':
        pathlist = config.get('Reading and Processing', 'List of raw data paths').split(', ')
        # checking the path list for valid directories
        rejectedpaths = [f for f in pathlist if os.path.isdir(f) is False]
        if len(rejectedpaths) != 0:
            for f in rejectedpaths:
                pathlist.remove(f)  # Removing those paths that were not found
                print('{} is not a usable MD directory'.format(f))
        elif len(rejectedpaths) == len(pathlist):
            pathlist = None
            print('No list of raw data paths found... exiting')
            exit()
    else:
        pathlist = None
        print('No list of raw data paths found... exiting')
        exit()

    # (B) Datalabels
    if config.get('Reading and Processing', 'datalabels') != '':
        datalabels = config.get('Reading and Processing', 'datalabels').split(', ')
    else:
        datalabels = None  # This will tell the readdynos function to generate generic datalabels raw_1...

    # (C) Ligand atom renaming
    if config.has_option('Reading and Processing', 'rename ligand') and \
            config.get('Reading and Processing', 'rename ligand') == '1':
        if (config.has_option('Reading and Processing', 'ligand renaming dictionary path') and
                os.path.isfile(config.get('Reading and Processing', 'ligand renaming dictionary path')) is True):
            ligatomdictpath = config.get('Reading and Processing', 'ligand renaming dictionary path')
            ligrename = True
        else:
            ligrename = False
            ligatomdictpath = None
            print('No suitable dictionary provided for ligand atom renaming')
    else:
        ligrename = False
        ligatomdictpath = None
    if config.has_option('Reading and Processing', 'ligand atom group dictionary'):
        ligandatomgroupdict = ast.literal_eval(config.get('Reading and Processing', 'ligand atom group dictionary'))
    else:
        ligandatomgroupdict = None

    # (D) Residue atom renaming
    if config.has_option('Reading and Processing', 'rename residues') and \
            config.get('Reading and Processing', 'rename residues') == '1':
        if (config.has_option('Reading and Processing', 'residue renaming dictionary path') and
                os.path.isfile(config.get('Reading and Processing', 'residue renaming dictionary path')) is True):
            resatomdictpath = config.get('Reading and Processing', 'residue renaming dictionary path')
            resrename = True
        else:
            resrename = False
            resatomdictpath = None
            print('No suitable dictionary provided for residue atom renaming')
    else:
        resrename = False
        resatomdictpath = None

    # (E) Selection statements: key residues, unwanted interactions, minimum interaction frequency
    if config.has_option('Reading and Processing', 'key residue list'):
        # converting all numbers to integers
        keyreslist = [int(x) for x in config.get('Reading and Processing', 'key residue list').split(', ')]
    else:
        keyreslist = None
    if config.has_option('Reading and Processing', 'exclude interactions'):
        # Transforming the string into a list of lists
        interactexclude = [list(x.strip().split('.')) for x in
                           config.get('Reading and Processing', 'exclude interactions').split(', ')]
    else:
        interactexclude = None
    if config.has_option('Reading and Processing', 'interaction frequency cutoff'):
        interactcutoff = float(config.get('Reading and Processing', 'interaction frequency cutoff'))
    else:
        interactcutoff = None

    # (F) Identical atoms
    if config.has_option('Reading and Processing', 'identical ligand atoms list'):
        # Transforming the string into a list of lists
        identliglist = [list(x.split('.')) for x in config.get('Reading and Processing', 'identical ligand atoms list').split(', ')]
    else:
        identliglist = None
    if config.has_option('Reading and Processing', 'identical residue atoms list'):
        identreslist = config.get('Reading and Processing', 'identical residue atoms list')
    else:
        identreslist = None

    # (G) Processing modus
    if config.has_option('Reading and Processing', 'modus') and config.get('Reading and Processing', 'modus') != '':
        dmode = config.get('Reading and Processing', 'modus')
        print('Working in {} modus'.format(dmode))
    else:
        print('Working in bc modus')
        dmode = 'bc'

else:
    pathlist = ligrename =  resrename = datalabels = keyreslist = interactexclude = ligandatomgroupdict = \
        resatomdictpath = interactcutoff = identliglist = ligatomdictpath = dmode = None
    print('Config file misses "Reading" section.')
    exit()

if config.has_section('Analysis'):
    # (G) Analysis settings
    if config.get('Analysis', 'fe-som-data paths') != '':
        feinclude = True  # Adding the heme iron - SoM distance to the dataframe
        # Input should be tested - like raw data paths
        feinpathlist = config.get('Analysis', 'fe-som-data paths').split(', ')

        # checking the fein path list for valid files
        rejectedfepaths = [f for f in feinpathlist if os.path.isfile(f) is False]
        if len(rejectedfepaths) != 0:
            for f in rejectedfepaths:
                feinpathlist.remove(f)  # Removing those fein paths that were not found
                print('{} is not a usable FESOM file'.format(f))
        elif len(rejectedfepaths) == len(pathlist):
            feinpathlist = None
            print('No list of usable fein paths found... exiting')
            exit()
    else:
        print('No list of fein paths found.')
        feinpathlist = None
        feinclude = False
        if config.has_section('Analysis', 'fe-som-label') and config.get('Analysis', 'fe-som-label') != '':
            fesomlabel = config.get('Analysis', 'fe-som-label')
        else:
            fesomlabel = None
    if config.has_option('Analysis', 'simulation time'):
        simtime = int(config.get('Analysis', 'simulation time'))
    else:
        simtime = None
else:
    simtime = feinpathlist = fesomlabel = None
    feinclude = False
#############################################################################################################
### this section probably needs some reworking, especially concerning saving figure WITHOUT showning them ###
#############################################################################################################
if config.has_section('Plotting'):
    if config.has_option('Plotting', 'plot heatmap') and config.get('Plotting', 'plot heatmap') == '1':
        heatmapplotting = True
    else:
        heatmapplotting = False
    if (config.has_option('Plotting', 'plot clustered heatmap') and config.get('Plotting',
                                                                               'plot clustered heatmap') == '1'):
        clustermapplotting = True
    else:
        clustermapplotting = False
else:
    heatmapplotting = clustermapplotting = False

if config.has_section('Writing'):
    # Checking for output path
    if (config.has_option('Writing', 'write directory') and
            os.path.isdir(config.get('Writing', 'write directory')) is True):
        writepath = config.get('Writing', 'write directory')
        # Boolean to specify, whether the rawdataframe should be saved as a pickle.
        if (config.has_option('Writing', 'pickle raw dataframes') and config.get('Writing',
                                                                                 'pickle raw dataframes') == '1'):
            rawpickle = True
        else:
            rawpickle = False
        if config.has_option('Writing', 'pickle final dataframes') and config.get('Writing',
                                                                                  'pickle final dataframes') == '1':
            pickleout = True
        else:
            pickleout = False
        if config.has_option('Writing', 'save plots') and config.get('Writing', 'save plots') == '1':
            savefig = True
        else:
            savefig = False
    else:
        print('No valid write directory provided')
        writepath = None
        rawpickle = pickleout = savefig = False
else:
    writepath = None
    rawpickle = pickleout = savefig = False

### Stage (i) ###
# Reading and processing dynophores raw data data to new dataframe structure.
# Working with either distance- (dc) or bit-based (bc) dataframes
if dmode == 'bc':
    bcdfdict = process_rawdata(pathlist, datalabels, keyreslist, interactexclude, ligrename, ligatomdictpath, resrename,
                               ligandatomgroupdict, resatomdictpath, interactcutoff, identliglist, dmode='bc', rawpickle=False)
    dcdfdict = None

elif dmode == 'dc':
    dcdfdict = process_rawdata(pathlist, datalabels, keyreslist, interactexclude, ligrename, ligatomdictpath, resrename,
                               ligandatomgroupdict, resatomdictpath, interactcutoff, identliglist, dmode='dc',  rawpickle=False)
    bcdfdict = None

else:
    dcdfdict = None
    bcdfdict = None

### Stage (ii) ###
# Analysing processed dataframes
# (1) Optionally Reading in Fe-SoM-distance data
bcfedfdict = {}
dcfedfdict = {}

if dmode == 'bc' and bcdfdict is not None:
    # Reading in Fe-SoM-distance data to identify binding trends
    for n, k in enumerate(bcdfdict.keys()):
        datalabel = replicalabel = k
        bcdf = bcdfdict[k]
        if feinclude is True and feinpathlist is not None:
            feinpath = feinpathlist[n]
            bcfedf = addfesom(bcdf, feinpath, fesomlabel=None, dftype='bc')
            bcfedfdict[k] = bcfedf
            if pickleout is True:
                pickler(bcfedf, datalabel, writepath, 'bcfedf')

        else:
            feinpath = None
            bcfedf = None
            if pickleout is True:
                pickler(bcdf, datalabel, writepath, 'bcdf')

elif dmode == 'dc' and dcdfdict is not None:
    for n, k in enumerate(dcdfdict.keys()):
        datalabel = replicalabel = k
        dcdf = dcdfdict[k]
        if feinclude is True and feinpathlist is not None:
            feinpath = feinpathlist[n]
            dcfedf = addfesom(dcdf, feinpath, fesomlabel=None, dftype='dc')
            bcfedfdict[k] = dcfedf
            if pickleout is True:
                pickler(dcfedf, datalabel, writepath, 'dcfedf')
        else:
            feinpath = None
            dcfedf = None
            if pickleout is True:
                pickler(dcdf, datalabel, writepath, 'dcdf')


    # (12) Testing out search of unique binding modes
    # tdf = bcdf.drop_duplicates()
    # sns.clustermap(bcfedf.drop_duplicates().T, cmap='Blues')
    # sns.clustermap(bcfedf.drop_duplicates().T, cmap='Blues', col_cluster=False)

    # Clustering the interactions based on the 'Jaccard index' = Tanimoto score

    # ... same but without color bar  # sns.clustermap(bcfedf.T, cmap='Blues', col_cluster=False, metric='jaccard').cax.set_visible(False)