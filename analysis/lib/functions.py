"""This script contains all functions used to import, analyse and plot dynophores results"""


import os
import re
import pickle
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from string import ascii_uppercase

def getmatchingfiles(pattern1, pattern2, list):
    """
    From a list of files get a new sorted list of files matching both patterns 1 and 2
    :param pattern1:
    :param pattern2:
    :param list:
    :return:
    """
    return [file for file in sorted(list) if pattern1 in file and pattern2 in file]


def readdynos(pathlist, labellist='', rawpickle=False):
    """
    Reading in a list of dynophores rawdata and returning an dictionary of rawdataframes

    pathlist - list of dynophores 'raw_data' directories
    labellist='' - list of labels to be used for the simulations
    rawpickle=False - Write out each rawframe as a pickle
    """
    # Initializing dictionary of dataframes
    rawdfdict = {}

    for i, directory in enumerate(pathlist):
        # Defining the dictionary key of raw dataframe in the raw dataframe dictionary (rawdfdict)
        if labellist != '':
            rawdfdictkey = labellist[i]
        else:
            rawdfdictkey = 'raw_{}'.format(i + 1)

        totallist = []  # Initializing list to gather info from the dynophores raw data to construct the raw dataframe
        ### Reading in dynophores raw data ###
        for f in os.listdir(directory):
            frame = -1
            if "_envPartner_" in f:
                name = f.split("superFeature_", 1)[1].split("%", 1)[0]  # file name without jobname and '%<percentage>.txt'
                type = name.split('[')[0]
                residue = re.findall('envPartner_(.*?)\[', name)[0]
                found = re.findall('\[(.*?)\]', name)  # regular expression search for substrings between [...]
                ligatoms = found[0].split(',')
                ligatoms = [int(x) for x in ligatoms]  # Translating to integers
                ligatoms = sorted(ligatoms)  # Sorting
                ligatoms = ', '.join(str(x) for x in ligatoms)  # list of atoms to string separated by ', '
                resatoms = found[1].split(',')  # getting ligand atom numbers as list
                resatoms = [int(x) for x in resatoms]  # transforming number into integers for subsequent sorting
                resatoms = sorted(resatoms)
                resatoms = ', '.join(str(x) for x in resatoms)
                resid = int(residue.split('_')[1])
                resname = residue.split('_')[0]
                #chainid = residue.split('_')[2]
                r = np.loadtxt(directory + '/' + f, delimiter=' ')
                for line in r:
                    dist = float(list(line)[0])
                    booli = int(list(line)[1])
                    frame += 1
                    totallist.append([frame, resid, resname, resatoms, ligatoms, type, dist, booli])

        header = ['frame', 'resid', 'resname', 'resatoms', 'ligatoms', 'type', 'dist', 'booli']  # Providing the column names for the dataframe
        rawdf = pd.DataFrame(totallist, columns=header)  # Creating a dataframe from totallist
        rawdfdict[rawdfdictkey] = rawdf  # Adding raw dataframe to the respective dictionary
        rawdf.sort_values(['frame', 'resid'], ascending=[True, True], inplace=True)  # Sorting on frame and then on residue id

        if rawpickle is True:
            rawdf.to_pickle(directory + 'rawdf.p', compression='gzip')  # If you want to save the rawdf for later

    return rawdfdict


def excludeinteract(rawdataframe, interactexcludelist):
    """
    Exclude the following cases from the dynophores raw dataframe
    Case1: Specific interaction type, all residues - example: 113.AR
    Case2: Specific residue, all interaction types- example 113.all
    Case3: Specific residue, specific interaction type - example all.AR

    returns dataframe without unwanted/exlcuded interactions

    rawdataframe - dynophores raw dataframe
    interactexclude - list of residue.interaction pairs
    """
    for i in interactexcludelist:
        if not '' in i:
            r = i[0]  # residue
            t = i[1]  # type of interaction
            # Case1: Specific interaction type, all residues
            if r == 'all' and t != 'all':
                rawdataframe = rawdataframe[rawdataframe['type'] != t]
                print('Removed all entries with interaction type of "{}"'.format(t))
            # Case2: Specific residue, all interaction types
            elif r != 'all' and t == 'all':
                rawdataframe = rawdataframe[rawdataframe['resid'] != int(r)]
                print('Removed all entries with residue id of "{}"'.format(r))
            # Case3: Specific residue, specific interaction type
            else:
                rawdataframe = rawdataframe[(rawdataframe['resid'] != int(r)) | (rawdataframe['type'] != t)]
                print(
                    'Removed all entries with the combination: residue id of "{0}" and interaction type "{1}"'.format(r, t))
        else:
            print('Exclude interactions input "{}" was not understood and will be omitted.'.format(i))

        return rawdataframe


def getligatomdict(pdbtuplelist, ligname):
    """Generate ligand atom label translation dictionary from list of tuples ('pdb-path', 'replica labels')"""
    ligatomdict = {}
    for i in pdbtuplelist:
        f = i[0]  # pdbfile
        replicalabel = i[1]  # replica label
        with open(f, 'r') as p:
            for line in p:
                if line.startswith('ATOM') and ligname in line:
                    atomno = line.split()[1]
                    atomlabel = line.split()[2]
                    if not atomlabel.startswith('H'):
                        ligatomdict[(atomno, replicalabel)] = atomlabel  # add to ligatomlabel
    return ligatomdict  # has format key = (atomid, replica label); value = 'atomname'


def lignum2name(string, dictionary, replicalabel):
    """Translate ligand atom labels in dynophores rawdata from atom ids to atom names.
    Requires dictionary for translation!"""
    return ','.join([dictionary.get(('{}'.format(i.replace(',', '')), replicalabel)) for i in string.split()])


def ligrenamer2atomname(rawdataframe, datalabel, ligatomdict):
    """
    Translate ligand atom labels in dynophores rawdata from atom ids to atom names.
    Requires translation dictionary!!

    returns dataframe with renamed ligand atom

    """
    newlignamecolumns = rawdataframe['ligatoms'].transform(lambda x: lignum2name(x, ligatomdict, datalabel))
    rawdataframe = rawdataframe.assign(ligatoms=newlignamecolumns)

    return rawdataframe


def ligrenamer2group(rawdataframe, ligatomgroupdict):
    """
    Translate ligand atom names into atom groups using the renaming dictionary provided in the configfile

    returns dataframe without unwanted/exlcuded interactions

    rawdataframe - dynophores raw dataframe
    ligatomgroupdict -

    """
    if ligatomgroupdict is not None:
        newcolumnsdict = {}  # Constructing dictionary for renaming
        for c in rawdataframe['ligatoms']:
            searchstring = c.replace(',', ', ')  # Getting same formatting for searching dictionary
            newline = ligatomgroupdict.get(searchstring)
            if newline is not None:
                newcolumnsdict[c] = c.replace('{}'.format(c),
                                              '{}'.format(newline))  # replace lig atom labels with atom group label
            else:
                newcolumnsdict[c] = c  # AKA do not rename these...

        # Renaming atoms if dictionary value available in two step procedure
        newlignamecolumn = [newcolumnsdict[item] if newcolumnsdict[item] is not None else item for item in
                            rawdataframe['ligatoms']]  # list comprehension with if, else
        rawdataframe = rawdataframe.assign(ligatoms=newlignamecolumn)

    return rawdataframe


def getresatomdict(pdbtuplelist, keyresidue=None):
    """
    Generate residue atom label translation dictionary from list of tuples ('pdb-path', 'replica labels')
    and optional limit to one key residue
    """
    waternames = ['T3P', 'HOH', 'T4P', 'SPC']  # Could be extended..
    resatomdict = {}
    # Key residue was provided and will be used
    if keyresidue is not None:
        for i in pdbtuplelist:
            f = i[0]  # pdbfile
            replicalabel = i[1]  # replica label
            with open(f, 'r') as p:
                for line in p:
                    # Read only atom entries, which belong to the residue
                    if (line.startswith('ATOM') and keyresidue == line.split()[5] and
                            str(line.split()[3]) not in waternames):
                        atomno = line.split()[1]
                        atomlabel = line.split()[2]
                        residue = str(line.split()[3]) + str(line.split()[5])
                        resatomdict[(atomno, replicalabel)] = [atomlabel, residue]  # add to resgatomlabel
    # all atoms will be read into a dictionary
    else:
        for i in pdbtuplelist:
            f = i[0]  # pdbfile
            replicalabel = i[1]  # replica label
            with open(f, 'r') as p:
                for line in p:
                    # Read only atom entries, which belong to the residue
                    if line.startswith('ATOM') and str(line.split()[3]) not in waternames:
                        atomno = line.split()[1]
                        atomlabel = line.split()[2]
                        residue = str(line.split()[3]) + str(line.split()[5])
                        resatomdict[(atomno, replicalabel)] = [atomlabel, residue]  # add to resgatomlabel
    return resatomdict


def resrenamer(dataframe, datalabel, renamedict, dftype='raw'):
    """
    Function to rename residue atoms in pandas DataFrame of dynophores analysis data
    to residuename.<list,of,residue,atoms>

    Datalabel - refers to the replica name that is also used by 'readdynophores'
    renamedict - can be obtained with getresatomdict and respective pdbtuplelist
    dftype='raw' - Type of dataframe; 'raw' or 'piv'
    """

    # Dealing with dynophores raw dataframe
    if dftype == 'raw':
        rawdataframe = dataframe
        newresatomcolumn = []
        for i in rawdataframe['resatoms']:
            newline = None
            # If there is only one residue atom involved
            if len(i.split(', ')) == 1:
                atomname, residue = renamedict.get((i, datalabel))
                if atomname is not None:
                    newresatomcolumn.append(str(atomname))  # adding renamed entry
                else:
                    newresatomcolumn.append(str('NAN'))
            # If there is more than one residue atom involved
            else:
                atomtranstmp = []
                newline = None
                for atom in i.split(', '):
                    atomname, residue = renamedict.get((atom, datalabel))
                    if atomname is not None:
                        atomtranstmp.append(str(atomname))
                    else:
                        atomtranstmp.append(str('NAN'))
                newline = ','.join(atomtranstmp)
                newresatomcolumn.append(str(newline))
        rawdataframe = rawdataframe.assign(resatoms=newresatomcolumn)
        print('Renamed residue atom numbers to residue atom names')

        return rawdataframe

    # Dealing with dynophores pivoted and processed dataframe
    elif dftype == 'piv':
        oldheader = dataframe.columns
        newheader = []

        for i in oldheader:
            i = i.split(';')  # old split column label

            # If there is only one residue atom involved
            if len(i[0].split(', ')) == 1:
                atomname = residue = None
                atomname, residue = renamedict[(i[0], datalabel)]  # getting information from the renaming dictionary
                newname = residue + '.' + atomname
                if newname is not None:
                    m = '; '.join([newname, i[1], i[2]])  # new column label
                    newheader.append(m)
                else:
                    print('Cannot rename {}'.format(i))
            # If there is more than one residue atom involved
            else:
                atomtranstmp = []
                for atom in i[0].split(', '):
                    atomname, residue = renamedict[(atom, datalabel)]
                    atomtranstmp.append((atomname, residue))
                residue = atomtranstmp[0][1]
                newname = residue + '.' + ','.join(
                    [x[0] for x in atomtranstmp])  # adding atom names of residue after dot
                if newname is not None:
                    m = '; '.join([newname, i[1], i[2]])  # new column label
                    newheader.append(m)
                else:
                    print('Cannot rename {}'.format(i))

        if len(dataframe.columns.tolist()) == len(newheader):
            return dataframe.rename(columns=dict(zip(dataframe.columns.tolist(), newheader)))

        else:
            print('Cannot finalize renaming, since old header and new header do not have the same size')
            return

    else:
        print('Cannot initiate renaming, since dataframe type was not unterstood')
        return


def interactfreqcutoff(dataframe, cutoff):
    for interact in dataframe.columns:
        if dataframe[interact].mean() * 100 < cutoff:
            dataframe.drop(labels=interact, axis=1, inplace=True)

    newinteractcount = len(list(dataframe.columns))  # new list without discarded interactions

    print('Discarding interactions with frequency below {} %'.format(cutoff))
    print('Reduced combinations of interactions to {}.'.format(newinteractcount))
    return dataframe


def mergeidentinteract(dataframe, identliglist):
    # Iterating over list of list of identical ligand atoms

    for identatoms in identliglist:
        mergedict = {}  # Initializing dictionary to collect columnnames to replace
        newlabel = identatoms[2]  # New atom label to be used after merging and renaming

        # Adding column names to dictionary if identatom is involved
        for interact in dataframe.columns:
            if identatoms[0] == interact.split(';')[-2]:
                newinteract = interact.replace('{}'.format(identatoms[0]), '{}'.format(newlabel))
                try:
                    mergedict[newinteract].append(interact)
                except KeyError:
                    mergedict[newinteract] = []
                    mergedict[newinteract].append(interact)
            elif identatoms[1] == interact.split(';')[-2]:
                newinteract = interact.replace('{}'.format(identatoms[1]), '{}'.format(newlabel))
                try:
                    mergedict[newinteract].append(interact)
                except KeyError:
                    mergedict[newinteract] = []
                    mergedict[newinteract].append(interact)

        # Iterating over dictionary keys and values
        for newinteract, oldinteracts in mergedict.items():
            # Check if there are actually two columns to merge
            if len(oldinteracts) == 2:
                # Creating the merged data
                newinteractdata = []
                for i, j in zip(dataframe[oldinteracts[0]], dataframe[oldinteracts[1]]):
                    if i == 1 or j == 1:
                        newinteractdata.append(1)
                    else:
                        newinteractdata.append(0)
                dataframe.drop(columns=oldinteracts, inplace=True)  # Dropping single columns
                dataframe['{}'.format(newinteract)] = newinteractdata  # Adding merged columns

    print('Merging interactions with identical identifier/ column name into one')
    print('Reduced combinations of interactions to {}.'.format(len(dataframe.columns)))
    return dataframe


def addfesom(dataframe, feinpath, fesomlabel=None, dftype='bc'):
    """
    Load in fesom distance data
    :param dataframein:
    :param dataframeout:
    :param feinpath:
    :param fesomlabel:
    :param dftype:
    :return:
    """
    fein = pd.read_csv(feinpath, usecols=[2], header=None, names=['fesomdist'], sep=',')
    fein.set_index(dataframe.index, inplace=True)
    # CAVE: pandas.DataFrame.where() is different from numpy.where() !!!
    if dftype == 'bc':
        fein.where(fein >= 6.0, 1, inplace=True)
        fein.where(fein < 6.0, 0, inplace=True)
        fein = fein.astype(dtype='int64')
        dataframe = dataframe.assign(fesomdist=fein['fesomdist'])
    elif dftype == 'dc':
        fein = fein.astype(dtype='float64')
        dataframe = dataframe.assign(fesomdist=fein['fesomdist'])
    if fesomlabel is not None:
        dataframe.rename(index=str, columns={"fesomdist": "{}".format(str(fesomlabel))}, inplace=True)

    print('Loaded in fesomdist data to new dataframe copy')
    return dataframe


def cm2inch(value):
    """Conversion from Centimeters to Inches used for figure size"""
    return value/2.54

def process_rawdata(pathlist, datalabels, keyreslist, interactexclude, ligrename, ligatomdictpath, resrename, ligandatomgroupdict,
                    resatomdictpath, interactcutoff, identliglist, dmode='bc',  rawpickle=False):
    """
    Reading in rawdata of dynophores analysis and processing it in the desired way.
    :param pathlist:
    :param datalabels:
    :param keyreslist:
    :param interactexclude:
    :param ligrename:
    :param ligatomdictpath:
    :param resrename:
    :param ligandatomgroupdict:
    :param resatomdictpath:
    :param interactcutoff:
    :param identliglist:
    :param dmode:
    :param rawpickle:
    :return:
    """

    bcdfdict = {}  # Initialising dictionary to collect the processed bit-based dataframes
    dcdfdict = {}  # Initialising dictionary to collect the processed distance-based dataframes

    # Reading in raw data into dictionary
    if datalabels is not None:
        rawdfdict = readdynos(pathlist, rawpickle=rawpickle, labellist=datalabels)
    else:
        rawdfdict = readdynos(pathlist, rawpickle=rawpickle)

    # Iterating over raw data frames in the dictionary
    for n, k in enumerate(rawdfdict.keys()):
        datalabel = replicalabel = k
        print('Starting to read {}'.format(k))
        print('Path: {}'.format(pathlist[n]))

        rawdf = None  # deleting old variable definition .. just to be sure

        # (1) Selecting only key residue interactions, if key residue list is not emtpy
        if keyreslist is not None:
            rawdf = rawdfdict[k][rawdfdict[k]['resid'].isin(keyreslist)]
        else:
            rawdf = rawdfdict[k]

        # (2) Excluding unwanted interactions as stated in the interactexclude list
        if interactexclude is not None:
            rawdf = excludeinteract(rawdf, interactexclude)

        # (3a) Translate ligand atom labels in dynophores rawdata from atom ids to atom names
        if ligrename is True and ligatomdictpath is not None:
            ligatomdict = pickle.load(open(ligatomdictpath, "rb"))
            rawdf = ligrenamer2atomname(rawdf, datalabel, ligatomdict)
            # (3b) Translate ligand atom names into atom groups using the renaming dictionary provided in the configfile
            rawdf = ligrenamer2group(rawdf, ligandatomgroupdict)

        # (4) Translate residue atom labels in dynophores rawdata from atom ids to atom names
        if resrename is True:
            keyatomdict = pickle.load(open(resatomdictpath, "rb"))
            rawdf = resrenamer(rawdf, datalabel, keyatomdict, dftype='raw')

        # (5) Adding new columns:
        # (5a) Column named 'interaction' with residue name, residue id, residue atoms, ligand atoms & interaction type
        rawdf = rawdf.assign(
            interaction=rawdf['resname'] + ';' + rawdf['resid'].astype(str) + ';' + rawdf['resatoms'] + ';' + rawdf[
                'ligatoms'] + ';' + rawdf['type'])
        # (5b) Column with the name of the system or replica
        rawdf = rawdf.assign(system=k)

        # Looking for the number of unique combinations of interactions
        uniqinteractnum = len(rawdf['interaction'].unique().tolist())
        print('Found {} unique combinations of interactions aka binding modes'.format(uniqinteractnum))

        # (6) Processing data to new data frame structure. Working with either distance- or bit-based data frame
        if dmode == 'bc':
            print('Starting to process {0} as {1}'.format(k, dmode))
            print('Path: {}'.format(pathlist[n]))

            # Creating a new dataframe format as barcode 'bcdf' for easier analysis
            tempbcdf = rawdf[['frame', 'interaction', 'booli']]
            #  Pivoting the dataframe aka 1D to 2D data
            bcdf = tempbcdf.pivot(index='frame', columns='interaction', values='booli')
            bcdf = bcdf.astype(int)  # setting items to integers

            # (7) Removing interactions occurring less frequent than the cutoff
            if interactcutoff is not None and 0 < interactcutoff < 100:
                bcdf = interactfreqcutoff(bcdf, interactcutoff)
            else:
                print('Provided interaction cutoff is out of range 0.0-100.0 %.')

            # (8) Merging interactions with identical atoms using list of identical atoms from config (identliglist)
            if identliglist != '' and identliglist != [] and identliglist is not None:

                bcdf = mergeidentinteract(bcdf, identliglist)

            # sorting columns alphabetically
            bcdf = bcdf.reindex(sorted(bcdf.columns), axis=1)
            # Adding processed dataframe to dictionary
            bcdfdict[k] = bcdf

        if dmode == 'dc':
            print('Starting to process {0} as {1}'.format(k, dmode))
            print('Path: {}'.format(pathlist[n]))

            # Creating a new dataframe format as distance-based 'dcdf' for easier analysis
            tempdcdf = rawdf[['frame', 'interaction', 'dist']]
            # Pivoting the dataframe aka 1D to 2D data
            dcdf = tempdcdf.pivot(index='frame', columns='interaction', values='dist')
            dcdfdict[k] = dcdf  # Adding processed dataframe to dictionary

    # output of function
    if dmode == 'bc':
        return bcdfdict
    elif dmode == 'dc':
        return dcdfdict

def pickler(dataframe, datalabel, writepath, dataframetype):
    """
    This function writes out the finalized dataframes as a pickle
    :param dataframe: dataframe to be pickled
    :param datalabel: datalabel to be used to name the file
    :param writepath: path to write to
    :param dataframetype: type of dataframe bc or dc + fe or not fe
    :return: nothing, but writes the pickled dataframe
    """

    pickledir = writepath.rstrip("/") + '/data-pickled'
    if not os.path.exists(pickledir):
        os.makedirs(pickledir)
    picklepath = "{0}/{1}-{2}.p".format(pickledir, datalabel, str(dataframetype))
    dataframe.to_pickle(picklepath, compression='gzip')
    print('Wrote out dataframe as pickle to {0} ({1} bytes)'.format(picklepath, os.path.getsize(picklepath)))

def readpickledf(pickleddfpath):
    """
    Reads in Pandas DataFrames that are stored as pickle with gzip compression.
    :param pickleddfpath:
    :return:
    """
    df = pd.read_pickle(pickleddfpath, compression='gzip')
    return df


def restrictinteract(dataframe, interactionlist, mode=''):
    oldcolumnslist = dataframe.columns.to_list()
    if interactionlist and interactionlist != []:
        newcolumnslist = []
        for c in oldcolumnslist:
            for i in interactionlist:
                # If only the interaction type should be restricted.
                if mode == 'type':
                    if c.endswith(';{}'.format(i)):
                        newcolumnslist.append(c)
                # If the full interaction (res + lig atoms + type) should be restricted.
                elif mode == 'full':
                    if c == i:
                        newcolumnslist.append(c)
        rmcolumns = [c for c in oldcolumnslist if c not in newcolumnslist]  # list of columns to be removed
        dataframe.drop(labels=rmcolumns, axis=1, inplace=True)
        return dataframe
    else:
        print('interaction list is empty or not usable.')


def groupatoms(dataframe, groupdict):
    oldheader = dataframe.columns.to_list()
    newheader = []
    for c in oldheader:
        a = c.split(';')[3]
        if a in groupdict.keys():
            newheader.append(c.replace('{}'.format(a), '{}'.format(groupdict[a])))
        else:
            newheader.append(c)
    return dataframe.rename(columns=dict(zip(dataframe.columns.tolist(), newheader)))


def dataframeanalysis(df, outcsvpath):
    """
    #### HEAVYLY UNDER CONSTRUCTION ###
    :param dataframe:
    :param outcsvpath:
    :return:
    """
    means = df.mean()  #

    # Following functionalities will be included!

    # bivariate analysis, see: /home/davidm/Tools/dynophores/dynos-analysis-readpickle-bivariate.py

    # Binding mode clustering?!

    # Plotting .. but put this in an separated function

    # get means and write these to CSV; see /home/davidm/Tools/dynophores/dynos-analysis-to-csv.py
    if outcsvpath:
        df.to_csv(outcsvpath, sep=',')  # might need some fixing concering header and columns, index and separator
