#!/usr/bin/env bash
##
## Transfers raw metagenomic data to data/ directory
##
## author: nlhuong90@gmail.com
## date: 12/20/2017

## Run the following commands from curie server where the data resides
## Specify username
USR=lanhuong
## Sepcify data directory to download
DIR=/scratch/users/lanhuong/Projects/PerturbationStudy/perturbation_16s/data/
## Would be nice if the data are gzipped...
rsync --compress --copy-links --ignore-existing -r /relman04/projects/hmd/MetaG/HMD_Mar2017/*Sub/*.fq $USR@dtn.sherlock.stanford.edu:$DIR/metagenomic
rsync --copy-links --ignore-existing -r /relman04/projects/hmd/MetaT/*Sub/*.fq.gz $USR@dtn.sherlock.stanford.edu:$DIR/metatranscriptomic
rsync --copy-links --ignore-existing -r /relman04/projects/hmd/MetaT/*Plate*/*.fq.gz $USR@dtn.sherlock.stanford.edu:$DIR/metatranscriptomic
