#!/usr/bin/env bash
##
## Pipline for processing metatranscriptomics data based on workflow:
## https://github.com/ParkinsonLab/2017-Microbiome-Workshop
##
## author: nlhuong90@gmail.com
## date: 2/18/2018

export n_threads=20
export ADD_BLAT=false
export INDEX_DB=true

## What Step to complete
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
export TAX_CLASS=false
export ASSEMBLE=false
export GENOME_ANN=false
export PROT_ANN=true


## Setup ------------

# Main directory
export STUDY_DIR=~/Projects/perturbation_16s
export MT_DIR=$STUDY_DIR/data/metatranscriptomics
# Reference directories
export REF_DIR=$STUDY_DIR/data/databases/
export KAIJUBD_DIR=$REF_DIR/kaijudb
# Input directories
export DATA_DIR=$MT_DIR/mouse
export INPUT_DIR=$DATA_DIR/input/
# Output directories
export OUTPUT_DIR=$DATA_DIR/output/
# Scripts and apps directories
export PYSCRIPT_DIR=$STUDY_DIR/metatranscriptomics/pyscripts_edited
export APP_DIR=~/.local/bin/

mkdir -p $DATA_DIR
mkdir -p $INPUT_DIR
mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR/QC

cd $OUTPUT_DIR

## Download data ------------

if $DOWNLOAD_DATA; then
    cd $INPUT_DIR
    wget https://github.com/ParkinsonLab/2017-Microbiome-Workshop/releases/download/Extra/precomputed_files.tar.gz
    tar -xvf precomputed_files.tar.gz mouse1.fastq
fi

## Check read quality -----------

if $QC; then
    cd $INPUT_DIR
    $APP_DIR/FastQC/fastqc mouse1.fastq
    mkdir -p $OUTPUT_DIR/QC
    mv *.html $OUTPUT_DIR/QC/
    mv *.zip $OUTPUT_DIR/QC/
    cd $OUTPUT_DIR
fi

## Remove adapter seqs and trim low-quality seqs  -----------

if $TRIM; then
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
    $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
       --fastq_filter mouse_merged_trim.fastq \
       --fastq_maxee 2.0 \
       --fastqout mouse_merged_qual.fastq
fi

if $QUAL_FLTR && ! $MERGE_PAIRS; then
    $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
       --fastq_filter mouse1_trim.fastq \
       --fastq_maxee 2.0 \
       --fastqout mouse1_qual.fastq
fi

## Remove duplicate reads ------------
if $RM_DUPL; then
    $APP_DIR/cdhit/cd-hit-auxtools/cd-hit-dup \
        -i mouse1_qual.fastq -o mouse1_unique.fastq
fi


## Remove vector contamination ------------
if $RM_VECTOR; then
    # generate an index for vector (UniVec_Core)  sequences
    if $INDEX_DB; then 
        bwa index -a bwtsw $REF_DIR/UniVec_Core
        samtools faidx $REF_DIR/UniVec_Core
        makeblastdb -in $REF_DIR/UniVec_Core -dbtype nucl
    fi
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
    
    $PYSCRIPT_DIR/1_BLAT_Filter.py \
        mouse1_univec_bwa.fastq \
        mouse1_univec.blatout \
        mouse1_univec_blat.fastq \
        mouse1_univec_blat_contaminats.fastq
fi

## Remove host reads ----------
# the same logic as the previous step
if $RM_HOST; then
    if $INDEX_DB; then
        bwa index -a bwtsw $REF_DIR/mouse_cds.fa
        samtools faidx $REF_DIR/mouse_cds.fa
        makeblastdb -in $REF_DIR/mouse_cds.fa -dbtype nucl
    fi
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
   
    $PYSCRIPT_DIR/1_BLAT_Filter.py \
        mouse1_mouse_bwa.fastq \
        mouse1_mouse.blatout \
        mouse1_mouse_blat.fastq \
        mouse1_mouse_blat_contaminats.fastq
fi

## Remove rRNA seqs -----------

if $RM_rRNA; then
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
    $PYSCRIPT_DIR/3_Reduplicate.py \
        mouse1_qual.fastq \
        mouse1_unique_mRNA.fastq \
        mouse1_unique.fastq.clstr \
        mouse1_mRNA.fastq
    $APP_DIR/FastQC/fastqc mouse1_mRNA.fastq
    mv *.html QC/
    mv *.zip QC/
fi

## Taxonomic Classification ------------

if $TAX_CLASS; then
   $APP_DIR/kaiju/bin/kaiju \
       -t $KAIJUBD_DIR/nodes.dmp \
       -f $KAIJUBD_DIR/kaiju_db.fmi \
       -i mouse1_mRNA.fastq \
       -z $n_threads \
       -o mouse1_classification.tsv

   $PYSCRIPT_DIR/4_Constrain_Classification.py \
       genus \
       mouse1_classification.tsv \
       $KAIJUBD_DIR/nodes.dmp \
       $KAIJUBD_DIR/names.dmp \
       mouse1_classification_genus.tsv

   $APP_DIR/kaiju/bin/kaijuReport \
       -t $KAIJUBD_DIR/nodes.dmp \
       -n $KAIJUBD_DIR/names.dmp \
       -i mouse1_classification_genus.tsv \
       -o mouse1_classification_Summary.txt \
       -r genus
fi

## Assemble reads into contigs -----------

if $ASSEMBLE; then
    # Build contigs
    $APP_DIR/SPAdes-3.11.1-Linux/bin/spades.py --rna \
        -s mouse1_mRNA.fastq \
        -o mouse1_spades
    mv mouse1_spades/transcripts.fasta mouse1_contigs.fasta
    
    # Use the contigs as the database and align unassembled mRNA reads to it
    bwa index -a bwtsw mouse1_contigs.fasta
    bwa mem -t $n_threads mouse1_contigs.fasta mouse1_mRNA.fastq > mouse1_contigs.sam
    
    # Make reads to contigs map 
    $PYSCRIPT_DIR/5_Contig_Map.py \
        mouse1_mRNA.fastq \
        mouse1_contigs.sam \
        mouse1_unassembled.fastq \
        mouse1_contigs_map.tsv
fi


## Genome annotation -----------

if $GENOME_ANN; then
    # Index DB
    if $INDEX_DB; then    
        bwa index -a bwtsw $REF_DIR/microbial_all_cds.fasta
        samtools faidx $REF_DIR/microbial_all_cds.fasta
    fi
    # Search DB
    bwa mem -t $n_threads $REF_DIR/microbial_all_cds.fasta \
        mouse1_contigs.fasta > mouse1_contigs_annotation_bwa.sam
    bwa mem -t $n_threads $REF_DIR/microbial_all_cds.fasta \
        mouse1_unassembled.fasta > mouse1_unassembled_annotation_bwa.sam
    
    if $ADD_BLAT; then
        # If we want to refine bwa alignments using bat
        # Extract reads that did not map to microbial genomes
        samtools view -bS \
           mouse1_contigs_annotation_bwa.sam > mouse1_contigs_annotation_bwa.bam
        samtools fastq -n -F 4 -0 \
            mouse1_contigs_annotation_aligned_bwa.fastq \
            mouse1_contigs_annotation_bwa.bam
        samtools fastq -n -f 4 -0 \
            mouse1_contigs_annotation_bwa.fastq \
            mouse1_contigs_annotation_bwa.bam
        $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
            --fastq_filter mouse1_contigs_annotation_bwa.fastq \
            --fastaout mouse1_contigs_annotation_bwa.fasta 
    
        # Additional alignments to reads not mapped by bwa
        blat $REF_DIR/microbial_all_cds.fasta \
            mouse1_contigs_annotation_bwa.fasta \
            -fine -q=rna -t=dna -out=blast8 \
            -noHead -minIdentity=90 -minScore=65  \
            mouse1_contigs_annotation.blatout

        $PYSCRIPT_DIR/1_BLAT_Filter.py \
            mouse1_contigs_annotation_bwa.fastq \
            mouse1_contigs_annotation.blatout \
            mouse1_contigs_annotation_blat.fastq \
            mouse1_contigs_annotation_blat_contaminats.fastq
    
        cat mouse1_contigs_annotation_aligned_bwa.fastq \
            mouse1_contigs_annotation_blat.fastq > \
            mouse1_contigs_annotation_bwa_blat.fastq
        ## DO THE SAME FOR UNASSEMBLED???
    fi

    $PYSCRIPT_DIR/6_BWA_Gene_Map.py \
        $REF_DIR/microbial_all_cds.fasta \
        mouse1_contigs_map.tsv \
        mouse1_genes_map.tsv \
        mouse1_genes.fasta \
        mouse1_contigs.fasta \
        mouse1_contigs_annotation_bwa.sam \
        mouse1_contigs_unmapped.fasta \
        mouse1_unassembled.fastq \
        mouse1_unassembled_annotation_bwa.sam \
        mouse1_unassembled_unmapped.fasta
fi
 
## Protein annotation -----------
if $PROT_ANN; then
    if $INDEX_DB; then 
       # Build index
       diamond makedb -p $n_threads --in $REF_DIR/nr -d $REF_DIR/nr
    fi
    mkdir -p dmnd_tmp
    diamond blastx -p $n_threads -d $REF_DIR/nr \
        -q mouse1_contigs_unmapped.fasta \
        -o mouse1_contigs.dmdout \
        -f 6 -t dmnd_tmp -k 10 \
        --id 85 --query-cover 65 --min-score 60
    diamond blastx -p $n_threads -d $REF_DIR/nr \
        -q mouse1_unassembled_unmapped.fasta \
        -o mouse1_unassembled.diamondout \
        -f 6 -t dmnd_tmp -k 10 \
        --id 85 --query-cover 65 --min-score 60
    
    $PYSCRIPT_DIR/7_Diamond_Protein_Map.py \
        $REF_DIR/nr \
        mouse1_contigs_map.tsv \
        mouse1_genes_map.tsv \
        mouse1_proteins.fasta \
        mouse1_contigs_unmapped.fasta \
        mouse1_contigs.dmdout \
        mouse1_contigs_unannotated.fasta \
        mouse1_unassembled_unmapped.fasta \
        mouse1_unassembled.dmdout \
        mouse1_unassembled_unannotated.fasta
fi


## Enzyme Annotation???
