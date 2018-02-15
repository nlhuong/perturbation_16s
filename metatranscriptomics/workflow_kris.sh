#!/usr/bin/env bash
##
## Metatranscriptomic processing of a single sample, following the workflow
## outlined in
## https://github.com/ParkinsonLab/2017-Microbiome-Workshop
##
## author: krissankaran@stanford.edu
## date: 2/14/2018

## download and organize the data / scripts
export WORK_DIR=$PWD
export MT_DIR=../data/metatranscriptomic/
export WF_DIR=$MT_DIR/workflow_raw/
export APP_DIR=~/.local/bin/
rm -rf $WF_DIR
mkdir $WF_DIR
mkdir $MT_DIR/QC/
wget https://github.com/ParkinsonLab/2017-Microbiome-Workshop/releases/download/Extra/precomputed_files.tar.gz

mv precomputed_files.tar.gz $WF_DIR
cd $WF_DIR
tar --wildcards -xvf precomputed_files.tar.gz *.py
tar -xvf precomputed_files.tar.gz mouse1.fastq

cd $WORK_DIR
mv $WF_DIR/*.py .

## generate the fastqc report
ln -fs $APP_DIR/Trimmomatic-0.36/adapters/TruSeq3-SE.fa Adapters
java -jar $APP_DIR/Trimmomatic-0.36/trimmomatic-0.36.jar \
     SE $WF_DIR/mouse1.fastq $WF_DIR/mouse1_trim.fastq \
     ILLUMINACLIP:Adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
$APP_DIR/FastQC/fastqc $WF_DIR/mouse1.fastq
$APP_DIR/FastQC/fastqc $WF_DIR/mouse1_trim.fastq
mv $WF_DIR/*.html $MT_DIR/QC/

## we'll have to run a paired end merging step in our analysis
## vsearch --fastq_mergepairs mouse1_trim.fastq --reverse mouse2_trim.fastq --fastqout mouse_merged_trim.fastq --fastqout_notmerged_fwd mouse1_merged_trim.fastq --fastqout_notmerged_rev mouse2_merged_trim.fastq

## global quality filtering
$APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
    --fastq_filter $WF_DIR/mouse1_trim.fastq --fastq_maxee 2.0 \
    --fastqout $WF_DIR/mouse1_qual.fastq

$APP_DIR/cdhit/cd-hit-auxtools/cd-hit-dup -i $WF_DIR/mouse1_qual.fastq -o $WF_DIR/mouse1_unique.fastq

## remove unwanted vector sequences
wget ftp://ftp.ncbi.nih.gov/pub/UniVec/UniVec_Core
mv UniVec_Core $WF_DIR
cd $WF_DIR

module load biology
module load bwa/0.7.17
module load samtools/1.6
module load ncbi-blast+/2.6.0
bwa index -a bwtsw UniVec_Core
samtools faidx UniVec_Core
makeblastdb -in UniVec_Core -dbtype nucl

bwa mem -t 4 UniVec_Core mouse1_unique.fastq > mouse1_univec_bwa.sam
samtools view -bS mouse1_univec_bwa.sam > mouse1_univec_bwa.bam
samtools fastq -n -F 4 -0 mouse1_univec_bwa_contaminats.fastq mouse1_univec_bwa.bam
samtools fastq -n -f 4 -0 mouse1_univec_bwa.fastq mouse1_univec_bwa.bam
$APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
    --fastq_filter mouse1_univec_bwa.fastq --fastaout mouse1_univec_bwa.fasta # hack to convert to fasta
blat -noHead -minIdentity=90 -minScore=65  UniVec_Core mouse1_univec_bwa.fasta -fine -q=rna -t=dna -out=blast8 mouse1_univec.blatout
$WORK_DIR/1_BLAT_Filter.py mouse1_univec_bwa.fastq mouse1_univec.blatout mouse1_univec_blat.fastq mouse1_univec_blat_contaminats.fastq
