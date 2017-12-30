#!/usr/bin/env bash
##
## Transfers raw metagenomic data to data/ directory
##
## author: sankaran.kris@gmail.com
## date: 11/27/2017

## Specify username
USR=lanhuong
## Sepcify data directory to download
DIR=../data/
## Would be nice if the data are gzipped...
rsync --compress --copy-links --ignore-existing -r $USR@curie.stanford.edu:/relman04/projects/hmd/MetaG/HMD_Mar2017/*Sub/*.fq $DIR/metagenomic
rsync --copy-links --ignore-existing -r $USR@curie.stanford.edu:/relman04/projects/hmd/MetaT/*Sub/*.fq.gz $DIR/metatranscriptomic
rsync --copy-links --ignore-existing -r $USR@curie.stanford.edu:/relman04/projects/hmd/MetaT/*Plate*/*.fq.gz $DIR/metatranscriptomic
