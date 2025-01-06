#!/bin/bash

# NOTE: This code was originally written by Gregory Matousek, and I adjusted some aspects of it so I could
# more easily submit my code to the farm to be analyzed. All credit in the slurm directory goes to him!!

VERSION="pass1"
TARGET_TYPES=(
    "ET"
    "NH3"
    "C"
    "CH2"
    "ND3"
    "CD2"
)
# Define sidisdvcs as an array of directory paths
sidisdvcs=(
    "/cache/clas12/rg-c/production/summer22/$VERSION/10.5gev/${TARGET_TYPES[0]}/dst/train/sidisdvcs/*"
    "/cache/clas12/rg-c/production/summer22/$VERSION/10.5gev/${TARGET_TYPES[1]}/dst/train/sidisdvcs/*"
    "/cache/clas12/rg-c/production/summer22/$VERSION/10.5gev/${TARGET_TYPES[2]}/dst/train/sidisdvcs/*"
    "/cache/clas12/rg-c/production/summer22/$VERSION/10.5gev/${TARGET_TYPES[3]}/dst/train/sidisdvcs/*"
    "/cache/clas12/rg-c/production/summer22/$VERSION/10.5gev/${TARGET_TYPES[4]}/dst/train/sidisdvcs/*"
    "/cache/clas12/rg-c/production/summer22/$VERSION/10.5gev/${TARGET_TYPES[5]}/dst/train/sidisdvcs/*"
)

destination=/w/hallb-scshelf2102/clas12/holmberg/Dilution_Factors/data_out
#rm -r $destination
#mkdir -p $destination

for i in {0..5}
do
    dirpath=${sidisdvcs[${i}]}
    for hipo in $dirpath
    do
	target=${TARGET_TYPES[${i}]}

        printf $hipo"    "$target 
	printf "\n"
        echo $(basename $hipo) | sed -e s/[^0-9]//g | while read -r line; 
        do
	    #printf ${dirpath} 
	    #printf "\n"    
            run="${line##0}"
	    echo $run
	    printf "\n"
	    sbatch /w/hallb-scshelf2102/clas12/holmberg/Dilution_Factors/slurm/DilutionSectorData.slurm $hipo $target $run
        done    
    done
done
