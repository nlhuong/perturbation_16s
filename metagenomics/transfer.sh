#!/usr/bin/env bash
##
## Transfers raw metagenomic data to data/ directory
##
## author: sankaran.kris@gmail.com
## date: 12/13/2017

## Specify username
USR=lanhuong
## Sepcify data directory to download
DIR=../data/
## Would be nice if the data are gzipped...
rsync --compress --copy-links --ignore-existing -r $USR@curie.stanford.edu:/relman04/projects/hmd/MetaG/HMD_Mar2017/*Sub/*.fq $DIR/metagenomic
<<<<<<< HEAD
rsync --copy-links --ignore-existing -r $USR@curie.stanford.edu:/relman04/projects/hmd/MetaT/*Sub/*.fq.gz $DIR/metatranscriptomic
rsync --copy-links --ignore-existing -r $USR@curie.stanford.edu:/relman04/projects/hmd/MetaT/*Plate*/*.fq.gz $DIR/metatranscriptomic
=======
rsync --compress --copy-links --ignore-existing -r $USR@curie.stanford.edu:/relman04/projects/hmd/MetaT/*Sub/*.fq.gz $DIR/metatranscriptomic
rsync --compress --copy-links --ignore-existing -r $USR@curie.stanford.edu:/relman04/projects/hmd/MetaT/*Plate/*.fq.gz $DIR/metatranscriptomic
>>>>>>> f290dbbc4d7df285c65ab628e6c693b04bcc8e3b
