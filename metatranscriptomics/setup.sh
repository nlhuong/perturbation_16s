#!/usr/bin/env bash
##
## Setup for metatranscriptomics data processing pipline.
## Mainly downloading reference databases.
## The pipeline is based on the workflow below:
## https://github.com/ParkinsonLab/2017-Microbiome-Workshop
##
## author: nlhuong90@gmail.com
## date: 2/19/2018

n_threads=10

## DIRECTORIES ----------
if [ -z ${SCRATCH+x} ]; then
    #the variable $SCRATCH is unset
    echo Working on ICME cluster
    module load gcc/gcc6
    BASE_DIR=~/Projects/perturbation_16s
    APP_DIR=~/.local/bin/
else
    echo Working on SHERLOCK cluster
    ## The following modules must be preloaded:
    module load python/2.7.13
    # module load py-biopython/1.70
    module load biology
    module load bwa/0.7.17
    module load samtools/1.6
    module load ncbi-blast+/2.6.0
    #module load python/3.6.1
    BASE_DIR=$SCRATCH/Projects/perturbation_16s
    APP_DIR=$SCRATCH/applications/bin/
fi

# Matetranscriptomics data folder
MT_DIR=$BASE_DIR/data/metatranscriptomics
# Pthon scripts
PYSCRIPT_DIR=$BASE_DIR/metatranscriptomics/pyscripts_edited
# Reference directories
REF_DIR=$PI_SCRATCH/resilience/databases
KAIJUBD_DIR=$REF_DIR/kaijudb
# SortMeRNA directories
SORTMERNA_DIR=$APP_DIR/sortmerna

mkdir -p $PYSCRIPT_DIR
mkdir -p $MT_DIR
mkdir -p $REF_DIR

cd $REF_DIR

## Download Scripts ------------
wget https://github.com/ParkinsonLab/2017-Microbiome-Workshop/releases/download/Extra/precomputed_files.tar.gz
tar --wildcards -xvf precomputed_files.tar.gz *.py
mv *.py $PYSCRIPT_DIR
rm precomputed_files.tar.gz

# Pyscripts needs editing if used with python3
# (A) 1_BLAT_Filter.py has inconsitent spaces/tabs
# (B) 2_Infernal_filter.py needs changing set() to dictionary {} as SeqRecord
#     type is no longer hashable in python3.
# (B) multiple files changing print to python3; need to change print statements
#
# Edit scripts and save in a new folder in the following location:
# $BASE_DIR/metatranscriptomics/pyscripts_edited


## Download reference databases ----------
# Vector sequences
cd $REF_DIR
wget ftp://ftp.ncbi.nih.gov/pub/UniVec/UniVec_Core

# Homo sapiense reference genome
wget ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz
gzip -d Homo_sapiens.GRCh38.cds.all.fa.gz
mv Homo_sapiens.GRCh38.cds.all.fa human_cds.fa

# Mouse reference genome
# wget ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz
# gzip -d Mus_musculus.GRCm38.cds.all.fa.gz
# mv Mus_musculus.GRCm38.cds.all.fa mouse_cds.fa


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

# Download Nonredundant (NR) protein database
wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
gunzip nr.gz

## The following is Adapted from SAMSA2

# Download NCBI RefSeq database:
echo -e "NOTE: The databases are up to 28GB and may require hours to download. Users may want to consider running this download overnight.\n"
echo "NOW DOWNLOADING NCBI REFSEQ DATABASE AT: "; date
wget --no-check-certificate "https://bioshare.bioinformatics.ucdavis.edu/bioshare/download/2c8s521xj9907hn/RefSeq_bac.fa"

# Download SEED Subsystems database:
echo "NOW DOWNLOADING SEED SUBSYSTEMS DATABASE AT: "; date
wget --no-check-certificate "https://bioshare.bioinformatics.ucdavis.edu/bioshare/download/2c8s521xj9907hn/subsys_db.fa" # Download non-redundant (NR) protein DB

# Download UniRef100 database:
wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz
gunzip uniref100.fasta.gz

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

# Index databases for DIAMOND
$APP_DIR/diamond makedb -p $n_threads --in $REF_DIR/RefSeq_bac.fa --db $REF_DIR/RefSeq_bac
$APP_DIR/diamond makedb -p $n_threads --in $REF_DIR/subsys_db.fa --db $REF_DIR/subsys_db
$APP_DIR/diamond makedb -p $n_threads --in $REF_DIR/nr -d $REF_DIR/nr
$APP_DIR/diamond makedb -p $n_threads --in $REF_DIR/uniref100.fasta -d $REF_DIR/uniref100

chmod g+xwr $REF_DIR/*

# Generated a text file with gene info [id, length, function/organism]
python $PYSCRIPT_DIR/db_to_pydict.py \
     $REF_DIR/RefSeq_bac.fa $REF_DIR/RefSeq_bac.tsv "refseq"
python $PYSCRIPT_DIR/db_to_pydict.py \
     $REF_DIR/nr $REF_DIR/nr.tsv "nr"
python $PYSCRIPT_DIR/db_to_pydict.py \
     $REF_DIR/subsys_db.fa $REF_DIR/subsys_db.tsv "seed"
python $PYSCRIPT_DIR/db_to_pydict.py \
     $REF_DIR/uniref100.fasta $REF_DIR/uniref100.tsv "uniref"


echo "Setup Completed!"
exit
