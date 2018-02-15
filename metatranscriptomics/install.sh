#!/usr/bin/env bash
##
## Installing things for the workflow in
## https://github.com/ParkinsonLab/2017-Microbiome-Workshop
##
## author: krissankaran@stanford.edu
## date: 2/14/2018

export APP_DIR=~/.local/bin/
export DB_DIR=/scratch/users/kriss1/research/perturbation_16s/data/kaijudb
cd $APP_DIR

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

## cd-hit
git clone git@github.com:weizhongli/cdhit.git
cd cdhit
make
cd cd-hit-auxtools
make

## BLAT
cd ../../
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat
chmod +x blat

## Infernal
wget http://eddylab.org/infernal/infernal-1.1.2-linux-intel-gcc.tar.gz
tar -zxvf infernal-1.1.2-linux-intel-gcc.tar.gz

## Kaiju
git clone https://github.com/bioinformatics-centre/kaiju.git
cd kaiju/src
make

cd ../..
mkdir $DB_DIR
cd $DB_DIR
$APP_DIR/kaiju/bin/makeDB.sh -r # make NCBI reference DB
cd ../../

## spades assembler
wget http://cab.spbu.ru/files/release3.11.1/SPAdes-3.11.1-Linux.tar.gz
tar -zxvf SPAdes-3.11.1-Linux.tar.gz
