#!/usr/bin/env bash
##
## Setup for metatranscriptomics data processing pipline. 
## Mainly downloading reference databases.
## The pipeline is based on the workflow below:
## https://github.com/ParkinsonLab/2017-Microbiome-Workshop
##
## author: nlhuong90@gmail.com
## date: 2/19/2018

export n_threads=10
# Location of directories
export STUDY_DIR=$SCRATCH/Projects/perturbation_16s
export PYSCRIPT_DIR=$STUDY_DIR/metatranscriptomics/pyscripts
export MT_DIR=$STUDY_DIR/data/metatranscriptomics
export REF_DIR=$STUDY_DIR/data/databases
export APP_DIR=$SCRATCH/applications/bin

mkdir -p $PYSCRIPT_DIR
mkdir -p $MT_DIR
mkdir -p $REF_DIR

cd $REF_DIR
## Download Scripts ------------

wget https://github.com/ParkinsonLab/2017-Microbiome-Workshop/releases/download/Extra/precomputed_files.tar.gz
tar --wildcards -xvf precomputed_files.tar.gz *.py
mv *.py $PYSCRIPT_DIR
rm precomputed_files.tar.gz

# Pyscripts needs editing if used with python 3
# file 1_(...) has inconsitent spaces/tabs 
# multiple files changing print to python3 convention
# rRNA removal step (file 2_Infernal_filter.py) is incompatible
# with python 3 as SeqRecord is not hashable edited file uses 
# dictionary instead of set(). Edited scripts are saved in
# a new folder
export PYSCRIPT_DIR=$STUDY_DIR/metatranscriptomics/pyscripts_edited

## Setup reference etc ----------

# vector sequences
cd $REF_DIR
wget ftp://ftp.ncbi.nih.gov/pub/UniVec/UniVec_Core

# Mouse reference genome
#wget ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz
#gzip -d Mus_musculus.GRCm38.cds.all.fa.gz
#mv Mus_musculus.GRCm38.cds.all.fa mouse_cds.fa

# Homo sapiense reference genome
wget ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz
gzip -d Homo_sapiens.GRCh38.cds.all.fa.gz
mv Homo_sapiens.GRCh38.cds.all.fa human_cds.fa

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
cd $REF_DIR

## Adapted from SAMSA2

# Download NCBI RefSeq database:
echo -e "NOTE: The databases are up to 28GB and may require hours to download. Users may want to consider running this download overnight.\n"
echo "NOW DOWNLOADING NCBI REFSEQ DATABASE AT: "; date
wget --no-check-certificate "https://bioshare.bioinformatics.ucdavis.edu/bioshare/download/2c8s521xj9907hn/RefSeq_bac.fa" 

# Download SEED Subsystems database:
echo "NOW DOWNLOADING SEED SUBSYSTEMS DATABASE AT: "; date
wget --no-check-certificate "https://bioshare.bioinformatics.ucdavis.edu/bioshare/download/2c8s521xj9907hn/subsys_db.fa" # Download non-redundant (NR) protein DB 
wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
gunzip nr.gz

# Download UniProt DB
wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz
tar -zxvf uniref100.fasta.gz
rm uniref100.fasta.gz

diamond makedb -p $n_threads --in $REF_DIR/RefSeq_bac.fa --db $REF_DIR/RefSeq_bac
diamond makedb -p $n_threads --in $REF_DIR/subsys_db.fa --db $REF_DIR/subsys_db
diamond makedb -p $n_threads --in $REF_DIR/nr -d $REF_DIR/nr
diamond makedb -p $n_threads --in $REF_DIR/uniref100.fasta -d $REF_DIR/uniref100

## Build database indexes -----------

# Index Vector DB
bwa index -a bwtsw $REF_DIR/UniVec_Core
samtools faidx $REF_DIR/UniVec_Core
makeblastdb -in $REF_DIR/UniVec_Core -dbtype nucl

# Index Host DB
bwa index -a bwtsw $REF_DIR/human_cds.fa
samtools faidx $REF_DIR/human_cds.fa
makeblastdb -in $REF_DIR/human_cds.fa -dbtype nucl

# Index Microbes DB
bwa index -a bwtsw $REF_DIR/microbial_all_cds.fasta
samtools faidx $REF_DIR/microbial_all_cds.fasta
makeblastdb -in $REF_DIR/microbial_all_cds.fasta -dbtype nucl

# Index SortMeRNA
$APP_DIR/sortmerna/build/Release/src/indexdb/indexdb \
    --ref $APP_DIR/sortmerna/rRNA_databases/silva-bac-16s-id90.fasta,$APP_DIR/sortmerna/index/silva-bac-16s -v

echo "Completed!"
exit

