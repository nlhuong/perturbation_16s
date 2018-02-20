#!/usr/bin/env bash
##
## Setup for metatranscriptomics data processing pipline. 
## Mainly downloading reference databases.
## The pipeline is based on the workflow below:
## https://github.com/ParkinsonLab/2017-Microbiome-Workshop
##
## author: lanhuong90@gmail.com
## date: 2/19/2018

# Location of directories
export STUDY_DIR=~/Projects/PerturbationStudy/perturbation_16s
export PYSCRIPT_DIR=$STUDY_DIR/metatranscriptomics/pyscripts
export MT_DIR=$STUDY_DIR/data/metatranscriptomics
export REF_DIR=$REF_DIR/references

mkdir -p $PYSCRIPT_DIR
mkdir -p $MT_DIR
mkdir -p $REF_DSTUDY

## Setup reference etc ----------

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
rm Rfam.cm.gz

# Genome and protein annotation dbs
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/all.ffn.tar.gz
mkdir microbial_references
tar -zxvf all.ffn.tar.gz -C microbial_references/
cd microbial_references
# merge genomes into one fasta
find . -name '*ffn' -exec cat {} \;> ../microbial_all_cds.fasta 
cd ../
rm -r microbial_references/
rm all.ffn.tar.gz

wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
gunzip nr.gz

# wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
# tar -zxvf uniref90.fasta.gz

## Download Scripts ------------

wget https://github.com/ParkinsonLab/2017-Microbiome-Workshop/releases/download/Extra/precomputed_files.tar.gz
tar --wildcards -xvf precomputed_files.tar.gz *.py
mv *.py $PYSCRIPT_DIR
rm precomputed_files.tar.gz

# Pyscripts needs editing
# file 1_(...) has inconsitent spaces/tabs 
# multiple files changing print to python3 convention
# rRNA removal step (file 2_Infernal_filter.py) is incompatible
# with python 3 as SeqRecord is not hashable edited file uses 
# dictionary instead of set(). Edited scripts are saved in
# a new folder
export PYSCRIPT_DIR=$STUDY_DIR/metatranscriptomics/pyscripts_edited

## Build database indexes -----------

# Index Vector DB
bwa index -a bwtsw $REF_DIR/UniVec_Core
samtools faidx $REF_DIR/UniVec_Core
makeblastdb -in $REF_DIR/UniVec_Core -dbtype nucl

# Index Host DB
bwa index -a bwtsw $REF_DIR/mouse_cds.fa
samtools faidx $REF_DIR/mouse_cds.fa
makeblastdb -in $REF_DIR/mouse_cds.fa -dbtype nucl

# Index Microbes DB
bwa index -a bwtsw $REF_DIR/microbial_all_cds.fasta
samtools faidx $REF_DIR/microbial_all_cds.fasta

# Index NonRedundant Proteins DB
diamond makedb -p $n_threads --in $REF_DIR/nr -d $REF_DIR/nr

