#!/usr/bin/env bash

# This script combines the bvec and bval files 
# from two DWI acquisitions into a single bvec and bval file.
# It's assumed that bvec and bval have the same name except for the extensions.
# Usage: combine_bvecs_and_bvals.sh DWI1.bvec DWI2.bvec

# argument checking
if [ $# -ne 2 ]; then
    echo "Usage: $0 DWI1.bvec DWI2.bvec"
    exit 1
fi
# check if the files exist
if [ ! -f "$1" ]; then
    echo "Error: File $1 not found!"
    exit 1
fi
if [ ! -f "$2" ]; then
    echo "Error: File $2 not found!"
    exit 1
fi


# file names without extensions
DWI1=${1%%.*}
DWI2=${2%%.*}

# combine the bvec files
cat $DWI1.bvec | tr -s ' ' > DWI1.clean
cat $DWI2.bvec | tr -s ' ' > DWI2.clean
paste -d ' ' DWI1.clean DWI2.clean > combined_${DWI1}_${DWI2}.bvec

# combine the bval files
cat $DWI1.bval | tr -s ' ' > DWI1.clean1
cat $DWI2.bval | tr -s ' ' > DWI2.clean1
paste -d ' ' DWI1.clean1 DWI2.clean1 > combined_${DWI1}_${DWI2}.bval

rm -rf *clean*
