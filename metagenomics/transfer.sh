#!/usr/bin/env bash
##
## Transfers raw metagenomic data to data/ directory
##
## author: sankaran.kris@gmail.com
## date: 11/27/2017

## Would be nice if the data are gzipped...
rsync --compress --copy-links --ignore-existing -r kriss1@curie.stanford.edu:/relman04/projects/hmd/MetaG/HMD_Mar2017/*Sub/*.fq ../data/
