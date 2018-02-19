#!/usr/bin/env bash
##
## Installing things for the workflow in
## https://github.com/ParkinsonLab/2017-Microbiome-Workshop
##
## author: nlhuong90@gmail.com
## date: 2/18/2018

export STUDY_DIR=~/Projects/PerturbationStudy/perturbation_16s
export APP_DIR=~/.local/bin/
cd $APP_DIR

# fastqc
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
unzip fastqc_v0.11.7.zip
rm fastqc_v0.11.7.zip
cd FastQC
chmod +x fastqc
echo  "export PATH=$APP_DIR/FastQC:\$PATH" >> ~/.bashrc

# trimmomatic
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
unzip Trimmomatic-0.36.zip
rm Trimmomatic-0.36.zip
chmod -R +x Trimmomatic-0.36
echo  "export PATH=$APP_DIR/Trimmomatic-0.36:\$PATH" >> ~/.bashrc

# vsearch
wget https://github.com/torognes/vsearch/releases/download/v2.7.0/vsearch-2.7.0-linux-x86_64.tar.gz
tar xzf vsearch-2.7.0-linux-x86_64.tar.gz
rm vsearch-2.7.0-linux-x86_64.tar.gz
chmod -R +x vsearch-2.7.0-linux-x86_64
echo  "export PATH=$APP_DIR/vsearch-2.7.0-linux-x86_64/bin:\$PATH" >> ~/.bashrc

# cd-hit
git clone https://github.com/weizhongli/cdhit.git
cd cdhit
make
cd cd-hit-auxtools
make
cd ../..

# BLAT
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat
chmod +x blat

# bwa
conda install -c bioconda bwa 

# samtools
conda install -c bioconda samtools

# blast
conda install -c bioconda blast

# Biopython for scripts
#conda install -c anaconda biopython #hashing SeqRecord issues
conda install -c conda-forge biopython #1.70
# ATTENTION
# The python scripts needed editing due to inconsistent tabbing
# and also print finction not compatible with python3

## Infernal
wget http://eddylab.org/infernal/infernal-1.1.2-linux-intel-gcc.tar.gz
tar -zxvf infernal-1.1.2-linux-intel-gcc.tar.gz
rm infernal-1.1.2-linux-intel-gcc.tar.gz

## Kaiju
git clone https://github.com/bioinformatics-centre/kaiju.git
cd kaiju/src
make

export DB_DIR=$STUDY_DIR/data/kaijudb
mkdir -p $DB_DIR
cd $DB_DIR
$APP_DIR/kaiju/bin/makeDB.sh -r # make NCBI reference DB
$APP_DIR/kaiju/bin/makeDB.sh -r --noDL
## Make custom library?
#$APP_DIR/kaiju/bin/mkbwt -o kaiju_db -nThreads 1 kaiju_db.faa
#$APP_DIR/kaiju/bin/mkfmi kaiju_db
cd $APP_DIR


## spades assembler
wget http://cab.spbu.ru/files/release3.11.1/SPAdes-3.11.1-Linux.tar.gz
tar -zxvf SPAdes-3.11.1-Linux.tar.gz
rm SPAdes-3.11.1-Linux.tar.gz

## Diamond
wget http://github.com/bbuchfink/diamond/releases/download/v0.9.17/diamond-linux64.tar.gz
tar -xzf diamond-linux64.tar.gz
rm diamond-linux64.tar.gz


source ~/.bashrc
