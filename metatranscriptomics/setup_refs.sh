#!/usr/bin/env bash
##
## Download all the reference databases required for the workflow
## https://github.com/ParkinsonLab/2017-Microbiome-Workshop
##
## author: krissankaran@stanford.edu
## date: 2/15/2018

rm -rf $PROCESS_DIR
mkdir $PROCESS_DIR
mkdir $REF_DIR
mkdir $MT_DIR/QC/

wget https://github.com/ParkinsonLab/2017-Microbiome-Workshop/releases/download/Extra/precomputed_files.tar.gz

mv precomputed_files.tar.gz $PROCESS_DIR
cd $PROCESS_DIR
tar --wildcards -xvf precomputed_files.tar.gz *.py
tar -xvf precomputed_files.tar.gz mouse1.fastq
rm precomputed_files.tar.gz

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
rm Rfam.cm.gz

# Genome and protein annotation dbs
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/all.ffn.tar.gz
mkdir microbial_references
tar -zxvf all.ffn.tar.gz -C microbial_references/
cd microbial_references
find . -name '*ffn' -exec cat {} \;> ../microbial_all_cds.fasta # merge genomes into one fasta
cd ../
rm -r microbial_references/
rm all.ffn.tar.gz

wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
gunzip nr.gz
rm nr.gz
