#!/usr/bin/env bash
##
## Download all the reference databases required for the workflow
## https://github.com/ParkinsonLab/2017-Microbiome-Workshop
##
## author: lanhuong90@gmail.com
## date: 2/18/2018

## What Step to complete
export SETUP=false
export DOWNLOAD_DATA=false
export QC=false
export TRIM=false
export MERGE_PAIRS=false
export QUAL_FLTR=false
export RM_DUPL=false
export RM_VECTOR=false
export RM_HOST=false
export RM_rRNA=false
export REREPLICATION=false

export n_threads=4
## Setup ------------

# Main directory
export STUDY_DIR=~/Projects/PerturbationStudy/perturbation_16s

# Input directories
export MT_DIR=$STUDY_DIR/data/metatranscriptomics
export REF_DIR=$MT_DIR/references
export DATA_DIR=$MT_DIR/mouse
export INPUT_DIR=$DATA_DIR/input/

# Output directories
export OUTPUT_DIR=$DATA_DIR/output/

# Script and app directories
export PYSCRIPT_DIR=$STUDY_DIR/metatranscriptomics/pyscripts
export APP_DIR=~/.local/bin/

mkdir -p $MT_DIR
mkdir -p $REF_DIR
mkdir -p $DATA_DIR
mkdir -p $INPUT_DIR
mkdir -p $OUTPUT_DIR
mkdir -p $PYSCRIPT_DIR


## Setup reference etc ----------

if $SETUP; then
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

    # wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
    # tar -zxvf uniref90.fasta.gz
fi


## Download data ------------

if $DOWNLOAD_DATA; then
    cd $INPUT_DIR
    wget https://github.com/ParkinsonLab/2017-Microbiome-Workshop/releases/download/Extra/precomputed_files.tar.gz
    tar -xvf precomputed_files.tar.gz mouse1.fastq
    tar --wildcards -xvf precomputed_files.tar.gz *.py
    mkdir -p $PYSCRIPT_DIR
    mv *.py $PYSCRIPT_DIR
    rm precomputed_files.tar.gz
fi
# Pyscripts needs editing firl 1_(...) inconsitent spaces/tabs and also print not compatible toth python3
export PYSCRIPT_DIR=$STUDY_DIR/metatranscriptomics/pyscripts_edited

## Check read quality -----------

if $QC; then
    cd $INPUT_DIR
    $APP_DIR/FastQC/fastqc mouse1.fastq
    mkdir -p $OUTPUT_DIR/QC
    mv *.html $OUTPUT_DIR/QC/
    mv *.zip $OUTPUT_DIR/QC/
fi

## Remove adapter seqs and trim low-quality seqs  -----------

if $TRIM; then
   cd $OUTPUT_DIR
   mkdir -p $OUTPUT_DIR/QC
   ln -fs $APP_DIR/Trimmomatic-0.36/adapters/TruSeq3-SE.fa Adapters
   java -jar $APP_DIR/Trimmomatic-0.36/trimmomatic-0.36.jar \
       SE $INPUT_DIR/mouse1.fastq $OUTPUT_DIR/mouse1_trim.fastq \
       ILLUMINACLIP:Adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
   $APP_DIR/FastQC/fastqc $OUTPUT_DIR/mouse1_trim.fastq
   mv $OUTPUT_DIR/*.html $OUTPUT_DIR/QC/
   mv $OUTPUT_DIR/*.zip $OUTPUT_DIR/QC/
fi

## Merge paired ends ------------
if $MERGE_PAIRS; then
    cd $OUTPUT_DIR
    mkdir -p $OUTPUT_DIR/QC
    ln -fs $APP_DIR/Trimmomatic-0.36/adapters/TruSeq3-SE.fa Adapters
    java -jar $APP_DIR/Trimmomatic-0.36/trimmomatic-0.36.jar \
        SE $INPUT_DIR/mouse2.fastq $OUTPUT_DIR/mouse2_trim.fastq \
        ILLUMINACLIP:Adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
    $APP_DIR/FastQC/fastqc $OUTPUT_DIR/mouse2_trim.fastq

    $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
        --fastq_mergepairs mouse1_trim.fastq \
        --reverse mouse2_trim.fastq \
        --fastqout mouse_merged_trim.fastq \
        --fastqout_notmerged_fwd mouse1_merged_trim.fastq \
        --fastqout_notmerged_rev mouse2_merged_trim.fastq

    $APP_DIR/FastQC/fastqc $OUTPUT_DIR/mouse_merged_trim.fastq
    mv $OUTPUT_DIR/*.html $OUTPUT_DIR/QC/
    mv $OUTPUT_DIR/*.zip $OUTPUT_DIR/QC/
fi


## Quality filtering ------------
if $QUAL_FLTR && $MERGE_PAIRS; then
    cd $OUTPUT_DIR
    $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
       --fastq_filter mouse_merged_trim.fastq \
       --fastq_maxee 2.0 \
       --fastqout mouse_merged_qual.fastq
fi

if $QUAL_FLTR && ! $MERGE_PAIRS; then
    cd $OUTPUT_DIR
    $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
       --fastq_filter mouse1_trim.fastq \
       --fastq_maxee 2.0 \
       --fastqout mouse1_qual.fastq
fi

## Remove duplicate reads ------------
if $RM_DUPL; then
    cd $OUTPUT_DIR
    $APP_DIR/cdhit/cd-hit-auxtools/cd-hit-dup \
        -i mouse1_qual.fastq -o mouse1_unique.fastq
fi


## Remove vector contamination ------------
if $RM_VECTOR; then
    cd $OUTPUT_DIR
    # generate an index for vector (UniVec_Core)  sequences
    bwa index -a bwtsw $REF_DIR/UniVec_Core
    samtools faidx $REF_DIR/UniVec_Core
    makeblastdb -in $REF_DIR/UniVec_Core -dbtype nucl
    
    # align reads with vector db  and filter any reads that align to it
    bwa mem -t $n_threads $REF_DIR/UniVec_Core mouse1_unique.fastq > mouse1_univec_bwa.sam
    samtools view -bS mouse1_univec_bwa.sam > mouse1_univec_bwa.bam
    samtools fastq -n -F 4 -0 mouse1_univec_bwa_contaminats.fastq mouse1_univec_bwa.bam
    samtools fastq -n -f 4 -0 mouse1_univec_bwa.fastq mouse1_univec_bwa.bam
    
    # additional alignments for the reads with BLAT to filter out any remaining
    # reads that align to our vector contamination database but BLAT only accepts 
    # fasta files so we have to convert our reads from fastq to fasta w/ vsearch
    $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
        --fastq_filter mouse1_univec_bwa.fastq \
        --fastaout mouse1_univec_bwa.fasta 
    
    blat $REF_DIR/UniVec_Core \
        mouse1_univec_bwa.fasta \
        -noHead -minIdentity=90 -minScore=65 \
        -fine -q=rna -t=dna -out=blast8 \
        mouse1_univec.blatout
    
    if [[ -s mouse1_univec.blatout ]] ; then
        $PYSCRIPT_DIR/1_BLAT_Filter.py \
            mouse1_univec_bwa.fastq \
            mouse1_univec.blatout \
            mouse1_univec_blat.fastq \
            mouse1_univec_blat_contaminats.fastq
    else 
        echo "mouse1_univec.blatout is empty"
        cp mouse1_univec_bwa.fastq mouse1_univec_blat.fastq
    fi
fi

## Remove host reads ----------

if $RM_HOST; then
    cd $OUTPUT_DIR
    ## exact same logic as removing vectors
    bwa index -a bwtsw $REF_DIR/mouse_cds.fa
    samtools faidx $REF_DIR/mouse_cds.fa
    makeblastdb -in $REF_DIR/mouse_cds.fa -dbtype nucl

    bwa mem -t $n_threads $REF_DIR/mouse_cds.fa mouse1_univec_blat.fastq > mouse1_mouse_bwa.sam
    samtools view -bS mouse1_mouse_bwa.sam > mouse1_mouse_bwa.bam
    samtools fastq -n -F 4 -0 mouse1_mouse_bwa_contaminats.fastq mouse1_mouse_bwa.bam
    samtools fastq -n -f 4 -0 mouse1_mouse_bwa.fastq mouse1_mouse_bwa.bam
    
    $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
        --fastq_filter mouse1_mouse_bwa.fastq \
        --fastaout mouse1_mouse_bwa.fasta

    blat $REF_DIR/mouse_cds.fa \
        mouse1_mouse_bwa.fasta \
        -noHead -minIdentity=90 -minScore=65 \    
        -fine -q=rna -t=dna -out=blast8 \
        mouse1_mouse.blatout
   
    if [[ -s mouse1_mouse.blatout ]] ; then
         $PYSCRIPT_DIR/1_BLAT_Filter.py \
             mouse1_mouse_bwa.fastq \
             mouse1_mouse.blatout \
             mouse1_mouse_blat.fastq \
             mouse1_mouse_blat_contaminats.fastq
    else
        echo "mouse1_mouse.blatout is empty"i
        cp mouse1_mouse_bwa.fastq mouse1_mouse_blat.fastq
    fi
fi

## Remove rRNA seqs -----------

if $RM_rRNA; then
    cd $OUTPUT_DIR 
    $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
        --fastq_filter \
        mouse1_mouse_blat.fastq \
        --fastaout \
        mouse1_mouse_blat.fasta
    $APP_DIR/infernal-1.1.2-linux-intel-gcc/binaries/cmsearch \
        -o mouse1_rRNA.log \
        --tblout mouse1_rRNA.infernalout \
        --anytrunc --rfam \
        -E 0.001 \
        $REF_DIR/Rfam.cm \
        mouse1_mouse_blat.fasta
    $PYSCRIPT_DIR/2_Infernal_Filter.py \
        mouse1_mouse_blat.fastq \
        mouse1_rRNA.infernalout \
        mouse1_unique_mRNA.fastq \
        mouse1_unique_rRNA.fastq
fi

## Rereplication (add back the repeated reads) ----------

if $REREPLICATION; then
    cd $OUTPUT_DIR
    $PYSCRIPT_DIR/3_Reduplicate.py \
        mouse1_qual.fastq \
        mouse1_unique_mRNA.fastq \
        mouse1_unique.fastq.clstr \
        mouse1_mRNA.fastq
    $APP_DIR/FastQC/fastqc mouse1_mRNA.fastq
    mv *.html QC/
    mv *.zip QC/
fi





