#! /bin/bash
# This script is used to get the names of non-standard residues (ligand, cofactor).
# example bash /mdspace/davidm/tools/getresnames.sh <CMS/file/path>

#input file to be investigated
cmsfile="$@"
#residue names to be excluded
resnameexclude="ACE|ALA|ARG|ASN|ASP|CL|CYS|GLN|GLU|GLY|HIE|HIP|HIS|ILE|LEU|LYS|MET|NA|NMA|PHE|PRO|SER|SPC|T3P|THR|TRP|TYR|VAL"
# command yielding a list of non-standard residue names (hopefully also your ligand)
reslist=$( grep -Eo "[0-9]\.[0-9]{5}\ \".{4}\"" $cmsfile | grep -Eo "\".{4}\"" | sort | uniq | sed 's#"##g' | grep -vE $resnameexclude )

echo $reslist
