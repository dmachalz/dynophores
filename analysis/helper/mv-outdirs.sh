#!/bin/bash
# Moving data from dynophores calculation one level up and deleting the generic folder

olddirs=$( find "$( realpath $1 )" -name "DynophoreApp*" -type d )
for olddir in $olddirs; do
  cd $olddir
    mv * ../
  cd ../
  if [ -z "$( find "$olddir" -mindepth 1 )" ]; then
    rm -rf $olddir
  fi
done
