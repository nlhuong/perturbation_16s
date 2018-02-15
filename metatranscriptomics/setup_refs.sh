#!/usr/bin/env bash
##
## Download all the reference databases required for the workflow
## https://github.com/ParkinsonLab/2017-Microbiome-Workshop
##
## author: krissankaran@stanford.edu
## date: 2/15/2018

export SCRIPT_DIR=/scratch/users/kriss1/research/perturbation_16s/metatranscriptomics/
export MT_DIR=$SCRIPT_DIR/../data/metatranscriptomic/
export REF_DIR=$MT_DIR/references
export PROCESS_DIR=$MT_DIR/workflow_raw/
rm -rf $PROCESS_DIR
mkdir $PROCESS_DIR
mkdir $REF_DIR
mkdir $MT_DIR/QC/

wget https://github.com/ParkinsonLab/2017-Microbiome-Workshop/releases/download/Extra/precomputed_files.tar.gz

mv precomputed_files.tar.gz $PROCESS_DIR
cd $PROCESS_DIR
tar --wildcards -xvf precomputed_files.tar.gz *.py
tar -xvf precomputed_files.tar.gz mouse1.fastq

cd $SCRIPT_DIR
mv $PROCESS_DIR/*.py .

# vector sequences
cd $REF_DIR
wget ftp://ftp.ncbi.nih.gov/pub/UniVec/UniVec_Core

# Mouse reference genome
wget ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz
gzip -d Mus_musculus.GRCm38.cds.all.fa.gz
mv Mus_musculus.GRCm38.cds.all.fa mouse_cds.fa

# Protein families
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/all.ffn.tar.gz
wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
tar -zxvf all.ff.nn.tar.gz
gunzip nr.gz
