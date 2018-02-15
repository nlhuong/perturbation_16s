#!/usr/bin/env bash
##
## Installing things for the workflow in
## https://github.com/ParkinsonLab/2017-Microbiome-Workshop
##
## author: krissankaran@stanford.edu
## date: 2/14/2018

cd ~/.local/bin/

module load java/1.8.0_131
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
unzip fastqc_v0.11.7.zip
chmod +x fastqc
export PATH=$PATH:~/.local/bin/FastQC
