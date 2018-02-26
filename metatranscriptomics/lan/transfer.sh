#!/usr/bin/env bash
##
## Transfers raw metagenomic data to data/ directory
##
## author: nlhuong90@gmail.com
## date: 2/25/2018

## Specify username
USR=lanhuong
## Sepcify data directory to download to
PROJDIR=/home/lanhuong/Projects/perturbation_16s
DATADIR=$PROJDIR/data/metatranscriptomics/resilience
mkdir -p $DATADIR
## Would be nice if the data are gzipped...
rsync --copy-links --ignore-existing -r $USR@curie.stanford.edu:/relman04/projects/hmd/MetaT/Second_Pilot/*Sub $DATADIR
rsync --copy-links --ignore-existing -r $USR@curie.stanford.edu:/relman04/projects/hmd/MetaT/GORIG/Second_Pilot/RNA_Plate* $DATADIR
rsync --copy-links --ignore-existing -r $USR@curie.stanford.edu:/relman04/projects/hmd/MetaT/NoIntervention_5/RNA_Plate* $DATADIR

#rsync -avzm --stats --safe-links --ignore-existing --dry-run --human-readable $USR@curie.stanford.edu:/relman04/projects/hmd/MetaT/Second_Pilot/*Sub $DATADIR > $DATADIR/transfer.log
## From curie server:
#rsync --copy-links --ignore-existing -r lanhuong@icme-share.stanford.edu:/home/lanhuong/Projects/perturbation_16s/data/metatranscriptomics/resilience

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
