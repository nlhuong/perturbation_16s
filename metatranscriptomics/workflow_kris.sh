#!/usr/bin/env bash
##
## Metatranscriptomic processing of a single sample, following the workflow
## outlined in
## https://github.com/ParkinsonLab/2017-Microbiome-Workshop
##
## author: krissankaran@stanford.edu
## date: 2/14/2018

## download and organize the data / scripts
export WORK_DIR=$PWD
export MT_DIR=../data/metatranscriptomic/
export WF_DIR=$MT_DIR/workflow_raw/
rm -rf $WF_DIR
mkdir $WF_DIR
mkdir $MT_DIR/QC/
wget https://github.com/ParkinsonLab/2017-Microbiome-Workshop/releases/download/Extra/precomputed_files.tar.gz

mv precomputed_files.tar.gz $WF_DIR
cd $WF_DIR
tar --wildcards -xvf precomputed_files.tar.gz *.py
tar -xvf precomputed_files.tar.gz mouse1.fastq

cd $WORK_DIR
mv $WF_DIR/*.py .

## generate the fastqc report
fastqc $WF_DIR/mouse1.fastq
mv $WF_DIR/*.html $MT_DIR/QC/
