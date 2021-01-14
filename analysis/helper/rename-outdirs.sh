#!/usr/bin/env bash
#change to folder of interest.. as supplied in the starting command after this script
cd $@

for sys in $( ls -d *); do
  for oldpathname in $( find "$sys" -type d -name "D*" | sort); do
    oldname=$( basename $oldpathname)
    newname=$( basename $( find "$oldpathname" -maxdepth 1 -name "*_CloudCenteroidDistances.txt" ) | sed "s#_CloudCenteroidDistances.txt##;s#$sys##;s#_##" )
    newpathname=$( echo $oldpathname | sed "s#$oldname#$newname#")
    mv $oldpathname $newpathname
  done
done
