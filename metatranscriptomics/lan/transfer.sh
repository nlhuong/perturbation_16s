#!/usr/bin/env bash
##
## Transfers raw metagenomic data to data/ directory
##
## author: nlhuong90@gmail.com
## date: 2/25/2018

## Run the following form the curie server where the data resides.
# data directory 
DATADIR=/relman04/projects/hmd/MetaT/
DESTINATION=/scratch/users/$USER/Projects/perturbation_16s/data/metatranscriptomics/resilience/input/

cd $DATADIR/Second_Pilot/
for dir in *; do 
    if [[ -d $dir ]] && [[ $dir = *"_Sub"* ]]; then 
        echo Copying directory $dir to $DESTINATION/$dir
        for file in ./$dir/*.fq.gz; do
            rsync --copy-links --ignore-existing -r \
                  $file $USER@dtn.sherlock.stanford.edu:$DESTINATION/$dir/ &
            sleep 1
        done
    fi
done

cd $DATADIR/NoIntervention_5/
for dir in *; do 
    if [[ -d $dir ]] && [[ $dir = *"RNA_Plate"* ]]; then 
        echo Copying directory $dir to $DESTINATION/$dir 
        for file in ./$dir/*.fq.gz; do
            rsync --copy-links --ignore-existing -r \
                  $file $USER@dtn.sherlock.stanford.edu:$DESTINATION/$dir/ &
            sleep 1
        done
    fi
done

#[lanhuong@curie Second_Pilot]$ ls DBUr_Sub/ | wc -l
#173
#[lanhuong@curie Second_Pilot]$ ls DBVr_Sub/ | wc -l
#213
#[lanhuong@curie Second_Pilot]$ ls EAQr_Sub/ | wc -l
#151
#[lanhuong@curie Second_Pilot]$ ls EAYr_Sub/ | wc -l
#105
#[lanhuong@curie Second_Pilot]$ ls EBFr_Sub/ | wc -l
#180
#[lanhuong@curie NoIntervention_5]$ ls RNA_Plate_11_TrM31/ | wc -l
#192
#[lanhuong@curie NoIntervention_5]$ ls RNA_Plate_9_TrM31/ | wc -l
#192


