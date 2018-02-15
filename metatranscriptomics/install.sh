#!/usr/bin/env bash
##
## Installing things for the workflow in
## https://github.com/ParkinsonLab/2017-Microbiome-Workshop
##
## author: krissankaran@stanford.edu
## date: 2/14/2018

cd ~/.local/bin/

## fastqc
module load java/1.8.0_131
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
unzip fastqc_v0.11.7.zip
chmod +x fastqc

## trimmomatic
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
unzip Trimmomatic-0.36.zip
chmod -R +x Trimmomatic-0.36

## vsearch
wget https://github.com/torognes/vsearch/releases/download/v2.7.0/vsearch-2.7.0-linux-x86_64.tar.gz
tar xzf vsearch-2.7.0-linux-x86_64.tar.gz
