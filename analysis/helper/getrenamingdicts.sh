#! /bin/bash
# This is a bash wrapper script for getrenamingdicts.py

#This takes the provided directory
WORKDIR="$1"
pdbsuffix='-in-noh2o.pdb'
pdblist=$( find "$WORKDIR" -name "*$pdbsuffix" -type f )
# creating an array


echo "Looking for pdbs in $WORKDIR"

# Creating a datalabel array from pdb files
datalabels=()
for pdb in $pdblist; do
  datalabel=$( basename $pdb | sed "s#$pdbsuffix##" )
  datalabels+=("$datalabel")
done

echo $pdblist
echo ""
echo "${datalabels[@]}"

exit
# Config file to be used with python helper script getrenamingdicts.py
CFG=''
# Executing the python script
python3 /mdspace/davidm/tools/dynophores/analysis/helper/getrenamingdicts.py $CFG
