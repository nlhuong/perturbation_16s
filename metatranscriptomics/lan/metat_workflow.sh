#!/usr/bin/env bash
##
## Pipline for processing metatranscriptomics data based on workflow:
## https://github.com/ParkinsonLab/2017-Microbiome-Workshop
##
## author: nlhuong90@gmail.com
## date: 2/18/2018

#module load python/2.7.13
#module load py-biopython/1.70
# module load biology
# module load bwa/0.7.17
# module load samtools/1.6
# module load ncbi-blast+/2.6.0
module load gcc/gcc6

## Input
export input_base=${1:-DBUr_Sub}
export input_fwd=${2:-M3303_DBUsw_2r_TrM31_1P.fq.gz}
export input_rev=${3:-M3303_DBUsw_2r_TrM31_2P.fq.gz}
export base=${input_fwd%_1P.fq.gz}
export fwd=${input_fwd%.fq.gz}
export rev=${input_rev%.fq.gz}

## Parameters
export paired=true
export n_threads=${4:-10}
export use_sortmerna=true
export add_blat=false

## What Step to complete
export INDEX_DB=false

export TRIM=false
export MERGE_PAIRS=false
export QUAL_FLTR=false
export RM_DUPL=false
export RM_VECTOR=false
export RM_HOST=false
export RM_rRNA=false
export REREPLICATION=true
export TAX_CLASS=true
export DIAMOND_REFSEQ=true
export DIAMOND_SEED=true
export ASSEMBLE=true
export GENOME_ANN=true
export PROT_ANN=true

## Setup ------------

# Input directories
export BASE_DIR=~/Projects/perturbation_16s
export DATA_DIR=$BASE_DIR/data/metatranscriptomics/resilience
export INPUT_DIR=$DATA_DIR/input/$input_base
# Output directories
export OUTPUT_DIR=$DATA_DIR/output/$input_base
# Reference directories
export REF_DIR=$BASE_DIR/data/databases
export KAIJUBD_DIR=$REF_DIR/kaijudb
# Scripts and apps directories
export APP_DIR=~/.local/bin/
export SORTMERNA_DIR=$APP_DIR/sortmerna
export PYSCRIPT_DIR=$BASE_DIR/metatranscriptomics/pyscripts_edited

mkdir -p $DATA_DIR
mkdir -p $INPUT_DIR
mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR/QC

cd $OUTPUT_DIR

start0=`date +%s`
## Generate an index
if $INDEX_DB; then
    echo "====================================================================================="
    echo "Indexing databases ..."
    start=`date +%s`
    # Vector (UniVec_Core)  sequences
    bwa index -a bwtsw $REF_DIR/UniVec_Core
    samtools faidx $REF_DIR/UniVec_Core
    makeblastdb -in $REF_DIR/UniVec_Core -dbtype nucl
    # Human genome
    bwa index -a bwtsw $REF_DIR/human_cds.fa
    samtools faidx $REF_DIR/human_cds.fa
    makeblastdb -in $REF_DIR/human_cds.fa -dbtype nucl
    # All RefSeq microbiome genomes
    bwa index -a bwtsw $REF_DIR/microbial_all_cds.fasta
    samtools faidx $REF_DIR/microbial_all_cds.fasta
    # Build diamond compatible rRNA_databases
    diamond makedb -p $n_threads --in $REF_DIR/RefSeq_bac.fa --db $REF_DIR/RefSeq_bac
    diamond makedb -p $n_threads --in $REF_DIR/subsys_db.fa --db $REF_DIR/subsys_db
    diamond makedb -p $n_threads --in $REF_DIR/nr -d $REF_DIR/nr
    diamond makedb -p $n_threads --in $REF_DIR/uniref100.fasta -d $REF_DIR/uniref100
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Indexing runtime: $runtime min" >> $OUTPUT_DIR/${base}_time.log
fi


## Remove adapter seqs and trim low-quality seqs  ------------
if $TRIM && ! $paired; then
   echo "====================================================================================="
   echo "Trimming and removing adapters ..."
   start=`date +%s`
   ln -fs $APP_DIR/Trimmomatic-0.36/adapters/TruSeq3-SE.fa Adapters
   java -jar $APP_DIR/Trimmomatic-0.36/trimmomatic-0.36.jar \
       SE $INPUT_DIR/$input_fwd $OUTPUT_DIR/${fwd}_trim.fastq \
       ILLUMINACLIP:Adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
   # Check reads quality
   $APP_DIR/FastQC/fastqc $INPUT_DIR/$input_fwd
   $APP_DIR/FastQC/fastqc $OUTPUT_DIR/${fwd}_trim.fastq
   mv $OUTPUT_DIR/*.html $OUTPUT_DIR/QC/
   mv $OUTPUT_DIR/*.zip $OUTPUT_DIR/QC/
   end=`date +%s`
   runtime=$(((end-start)/60))
   echo "Trimming and QC runtime: $runtime min" >> $OUTPUT_DIR/${base}_time.log
fi


## Do the same for paired ends ------------
if $TRIM && $paired; then
    echo "====================================================================================="
    echo "Trimming and removing adapters ..."
    start=`date +%s`
    ln -fs $APP_DIR/Trimmomatic-0.36/adapters/TruSeq3-SE.fa Adapters
    java -jar $APP_DIR/Trimmomatic-0.36/trimmomatic-0.36.jar \
        PE $INPUT_DIR/$input_fwd $INPUT_DIR/$input_rev \
        $OUTPUT_DIR/${fwd}_paired_trim.fq $OUTPUT_DIR/${fwd}_unpaired_trim.fq \
        $OUTPUT_DIR/${rev}_paired_trim.fq $OUTPUT_DIR/${rev}_unpaired_trim.fq \
        ILLUMINACLIP:Adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
    # check reads quality
    $APP_DIR/FastQC/fastqc $INPUT_DIR/$input_fwd
    $APP_DIR/FastQC/fastqc $INPUT_DIR/$input_rev
    $APP_DIR/FastQC/fastqc $OUTPUT_DIR/${fwd}_paired_trim.fq
    $APP_DIR/FastQC/fastqc $OUTPUT_DIR/${rev}_paired_trim.fq
    mv $OUTPUT_DIR/*.html $OUTPUT_DIR/QC/
    mv $OUTPUT_DIR/*.zip $OUTPUT_DIR/QC/
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Trimming and QC runtime: $runtime min" >> $OUTPUT_DIR/${base}_time.log
fi


## Merge pairs ------------
if $MERGE_PAIRS; then
    echo "====================================================================================="
    echo "Merging paired reads ..."
    start=`date +%s`
    $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
        --fastq_mergepairs $OUTPUT_DIR/${fwd}_paired_trim.fq \
        --reverse $OUTPUT_DIR/${rev}_paired_trim.fq \
        --fastqout $OUTPUT_DIR/${base}_trim.fq \
        --fastqout_notmerged_fwd $OUTPUT_DIR/${base}_unmerged_trim_fwd.fq \
        --fastqout_notmerged_rev $OUTPUT_DIR/${base}_unmerged_trim_rev.fq \

    $APP_DIR/FastQC/fastqc $OUTPUT_DIR/${base}_trim.fq
    mv $OUTPUT_DIR/*.html $OUTPUT_DIR/QC/
    mv $OUTPUT_DIR/*.zip $OUTPUT_DIR/QC/
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Merging reads runtime: $runtime min" >> $OUTPUT_DIR/${base}_time.log
fi


## Global quality filtering ------------
if $QUAL_FLTR; then
    echo "====================================================================================="
    echo "Quality filtering ..."
    start=`date +%s`
    $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
       --fastq_filter $OUTPUT_DIR/${base}_trim.fq \
       --fastq_maxee 2.0 \
       --fastqout $OUTPUT_DIR/${base}_qual.fq
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Quality filtering runtime: $runtime min" >> $OUTPUT_DIR/${base}_time.log
fi


## Remove duplicate reads ------------
if $RM_DUPL; then
    echo "====================================================================================="
    echo "Remove duplicates ..."
    start=`date +%s`
    $APP_DIR/cdhit/cd-hit-auxtools/cd-hit-dup \
        -i $OUTPUT_DIR/${base}_qual.fq \
        -o $OUTPUT_DIR/${base}_unique.fq
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Remove duplicates runtime: $runtime min" >> $OUTPUT_DIR/${base}_time.log
fi


## Remove unwanted vector contamination ------------
if $RM_VECTOR; then
    echo "====================================================================================="
    echo "Remove vector sequences ..."
    start=`date +%s`
    # align reads with vector db  and filter any reads that align to it
    bwa mem -t $n_threads $REF_DIR/UniVec_Core \
        ${base}_unique.fq > ${base}_univec_bwa.sam
    samtools view -bS ${base}_univec_bwa.sam > ${base}_univec_bwa.bam
    samtools fastq -n -F 4 -0 ${base}_univec_bwa_contam.fq ${base}_univec_bwa.bam
    samtools fastq -n -f 4 -0 ${base}_univec_bwa.fq ${base}_univec_bwa.bam

    # additional alignments for the reads with BLAT to filter out any remaining
    # reads that align to our vector contamination database but BLAT only accepts
    # fasta files so we have to convert our reads from fastq to fasta w/ vsearch
    $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
        --fastq_filter ${base}_univec_bwa.fq \
        --fastaout ${base}_univec_bwa.fasta

    blat $REF_DIR/UniVec_Core \
        ${base}_univec_bwa.fasta \
        -noHead -minIdentity=90 -minScore=65 \
        -fine -q=rna -t=dna -out=blast8 \
        ${base}_univec.blatout

    $PYSCRIPT_DIR/1_BLAT_Filter.py \
        ${base}_univec_bwa.fq \
        ${base}_univec.blatout \
        ${base}_univec_blat.fq \
        ${base}_univec_blat_contam.fq

    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Removing vectors runtime: $runtime min" >> $OUTPUT_DIR/${base}_time.log
fi

## Remove host reads ----------
# the same logic as the previous step
if $RM_HOST; then
    echo "====================================================================================="
    echo "Remove host sequences ..."
    start=`date +%s`
    bwa mem -t $n_threads $REF_DIR/human_cds.fa \
        ${base}_univec_blat.fq > ${base}_human_bwa.sam
    samtools view -bS ${base}_human_bwa.sam > ${base}_human_bwa.bam
    samtools fastq -n -F 4 -0 ${base}_human_bwa_contam.fq ${base}_human_bwa.bam
    samtools fastq -n -f 4 -0 ${base}_human_bwa.fq ${base}_human_bwa.bam

    $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
        --fastq_filter ${base}_human_bwa.fq \
        --fastaout ${base}_human_bwa.fasta

    blat $REF_DIR/human_cds.fa \
        ${base}_human_bwa.fasta \
        -noHead -minIdentity=90 -minScore=65 \
        -fine -q=rna -t=dna -out=blast8 \
        ${base}_human.blatout

    $PYSCRIPT_DIR/1_BLAT_Filter.py \
        ${base}_human_bwa.fq \
        ${base}_human.blatout \
        ${base}_human_blat.fq \
        ${base}_human_blat_contam.fq

    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Removing host runtime: $runtime min" >> $OUTPUT_DIR/${base}_time.log
fi

## Remove rRNA seqs -----------
if $RM_rRNA && $use_sortmerna; then
    echo "====================================================================================="
    echo "Ribodepletion with SortMeRNA ..."
    start=`date +%s`
    $SORTMERNA_DIR/build/Release/src/sortmerna/sortmerna \
        --ref $SORTMERNA_DIR/rRNA_databases/silva-bac-16s-id90.fasta,$SORTMERNA_DIR/index/silva-bac-16s \
        --reads ${base}_human_blat.fq \
        --aligned ${base}_unique_rRNA \
        --other ${base}_unique_mRNA \
        --fastx --num_alignments 0 --log \
        -a $n_threads -v 
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "SortMeRNA ribodepletion runtime: $runtime min" >> $OUTPUT_DIR/${base}_time.log
fi

if $RM_rRNA && ! $use_sortmerna; then
    echo "====================================================================================="
    echo "Ribodepletion with infernalout ..."
    start=`date +%s`
    
    $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
        --fastq_filter ${base}_human_blat.fq \
        --fastaout ${base}_human_blat.fasta

    $APP_DIR/infernal-1.1.2-linux-intel-gcc/binaries/cmsearch \
        -o ${base}_rRNA.log \
        --tblout ${base}_rRNA.infernalout \
        --anytrunc --rfam -E 0.001 --cpu $n_threads \
        $REF_DIR/Rfam.cm \
        ${base}_human_blat.fasta
    
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Infernal ribodepletion runtime: $runtime min" >> $OUTPUT_DIR/${base}_time.log

    start=`date +%s`
    $PYSCRIPT_DIR/2_Infernal_Filter.py \
        ${base}_human_blat.fq \
        ${base}_rRNA.infernalout \
        ${base}_unique_mRNA.fq \
        ${base}_unique_rRNA.fq
    
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Python infernal filter runtime: $runtime min" >> $OUTPUT_DIR/${base}_time.log
fi


## Rereplication (add back the repeated reads) ----------
if $REREPLICATION; then
    echo "====================================================================================="
    echo "Rereplication ..."
    start=`date +%s`
    $PYSCRIPT_DIR/3_Reduplicate.py \
        ${base}_qual.fq \
        ${base}_unique_mRNA.fq \
        ${base}_unique.fq.clstr \
        ${base}_mRNA.fq

    $APP_DIR/FastQC/fastqc ${base}_mRNA.fq
    mv $OUTPUT_DIR/*.html $OUTPUT_DIR/QC/
    mv $OUTPUT_DIR/*.zip $OUTPUT_DIR/QC/
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Rereplication runtime: $runtime min" >> $OUTPUT_DIR/${base}_time.log
fi

## Taxonomic Classification ------------
if $TAX_CLASS; then
    echo "====================================================================================="
    echo "Taxonomy classification ..."
    start=`date +%s`
    $APP_DIR/kaiju/bin/kaiju \
        -t $KAIJUBD_DIR/nodes.dmp \
        -f $KAIJUBD_DIR/kaiju_db.fmi \
        -i ${base}_mRNA.fq \
        -z $n_threads \
        -o ${base}_tax_class.tsv
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Kaiju runtime: $runtime min" >> $OUTPUT_DIR/${base}_time.log
    
    start=`date +%s`
    $PYSCRIPT_DIR/4_Constrain_Classification.py \
        genus \
        ${base}_tax_class.tsv \
        $KAIJUBD_DIR/nodes.dmp \
        $KAIJUBD_DIR/names.dmp \
        ${base}_genus_class.tsv
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Python script process taxonomy result runtime: $runtime min" >> $OUTPUT_DIR/${base}_time.log

    start=`date +%s`
    $APP_DIR/kaiju/bin/kaijuReport \
        -t $KAIJUBD_DIR/nodes.dmp \
        -n $KAIJUBD_DIR/names.dmp \
        -i ${base}_genus_class.tsv \
        -o ${base}_genus_class_summary.txt \
        -r genus
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Aggregate taxonomy res runtime: $runtime min" >> $OUTPUT_DIR/${base}_time.log
fi


## DIAMOND annotation -----------
if $DIAMOND_REFSEQ; then
    echo "====================================================================================="
    echo "Refseq annotation with DIAMOND ..."
    start=`date +%s`
    # Align with RefSeq database
    mkdir -p dmnd_tmp
    diamond blastx -p $n_threads -db $REF_DIR/RefSeq_bac \
        -q ${base}_mRNA.fq -a ${base}.RefSeq \
        -t ./dmnd_tmp -k 1 --sensitive
    diamond view --daa ${base}.RefSeq.daa -f 6 \
        -o ${base}_refseq.dmdout

    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Refseq annotation (diamond) runtime: $runtime min" >> \
            $OUTPUT_DIR/${base}_time.log
fi


if $DIAMOND_SEED; then
    echo "====================================================================================="
    echo "SEED annotation with DIAMOND ..."
    start=`date +%s`
    # Align with SEED subsystem database
    mkdir -p dmnd_tmp
    diamond blastx -p $n_threads -db $REF_DIR/subsys_d \
        -q ${base}_mRNA.fq -a ${base}.Subsys \
        -t ./dmnd_tmp -k 1 --sensitive
    diamond view --daa ${base}.Subsys.daa -f 6 \
        -o ${base}_seed.dmdout
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "SEED annotation (diamond) runtime: $runtime min" >> \
        $OUTPUT_DIR/${base}_time.log
fi


## Assemble reads into contigs -----------
if $ASSEMBLE; then
    echo "====================================================================================="
    echo "Contigs assembly ..."
    start=`date +%s`
    # Build contigs
    $APP_DIR/SPAdes-3.11.1-Linux/bin/spades.py \
        --rna -t $n_threads \
        -s ${base}_mRNA.fq \
        -o ${base}_spades
    mv ${base}_spades/transcripts.fasta ${base}_contigs.fasta

    # Use the contigs as the database and align unassembled mRNA reads to it
    bwa index -a bwtsw ${base}_contigs.fasta
    bwa mem -t $n_threads ${base}_contigs.fasta ${base}_mRNA.fq > ${base}_contigs.sam

    # Make reads to contigs map
    $PYSCRIPT_DIR/5_Contig_Map.py \
        ${base}_mRNA.fq \
        ${base}_contigs.sam \
        ${base}_unassembled.fq \
        ${base}_contigs_map.tsv
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Assemble runtime: $runtime min" >> $OUTPUT_DIR/${base}_time.log
fi


## Genome annotation -----------
if $GENOME_ANN; then
    echo "====================================================================================="
    echo "Genome annotation with BWA ..."
    start=`date +%s`
    # BWA utilizes nucleotide searches, we rely on a microbial genome database
    # Search DB
    bwa mem -t $n_threads $REF_DIR/microbial_all_cds.fasta \
        ${base}_contigs.fasta > ${base}_contigs_annotation_bwa.sam
    bwa mem -t $n_threads $REF_DIR/microbial_all_cds.fasta \
        ${base}_unassembled.fq > ${base}_unassembled_annotation_bwa.sam

    # IS THE FOLLOWING EVEN USEFUL? HOW TO ADD TO 6_BWA_Gene_Map.py ???
    if $add_blat; then
        # If we want to refine bwa alignments using blat
        # Extract reads that did not map to microbial genomes
        samtools view -bS \
           ${base}_contigs_annotation_bwa.sam > ${base}_contigs_annotation_bwa.bam
        samtools fastq -n -F 4 -0 \
            ${base}_contigs_annotation_aligned_bwa.fq \
            ${base}_contigs_annotation_bwa.bam
        samtools fastq -n -f 4 -0 \
            ${base}_contigs_annotation_bwa.fq \
            ${base}_contigs_annotation_bwa.bam

        $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
            --fastq_filter ${base}_contigs_annotation_bwa.fq \
            --fastaout ${base}_contigs_annotation_bwa.fasta

        # Additional alignments to reads not mapped by bwa
        blat $REF_DIR/microbial_all_cds.fasta \
            ${base}_contigs_annotation_bwa.fasta \
            -fine -q=rna -t=dna -out=blast8 \
            -noHead -minIdentity=90 -minScore=65  \
            ${base}_contigs_annotation_bwa.blatout

        $PYSCRIPT_DIR/1_BLAT_Filter.py \
            ${base}_contigs_annotation_bwa.fq \
            ${base}_contigs_annotation_bwa.blatout \
            ${base}_contigs_annotation_blat.fq \
            ${base}_contigs_annotation_aligned_blat.fq

        cat ${base}_contigs_annotation_aligned_bwa.fq \
            ${base}_contigs_annotation_aligned_blat.fq > \
            ${base}_contigs_annotation_aligned_bwa_blat.fq
        ## DO THE SAME FOR UNASSEMBLED???
    fi

    $PYSCRIPT_DIR/6_BWA_Gene_Map.py \
        $REF_DIR/microbial_all_cds.fasta \
        ${base}_contigs_map.tsv \
        ${base}_genes_map.tsv \
        ${base}_genes.fasta \
        ${base}_contigs.fasta \
        ${base}_contigs_annotation_bwa.sam \
        ${base}_contigs_unmapped.fasta \
        ${base}_unassembled.fq \
        ${base}_unassembled_annotation_bwa.sam \
        ${base}_unassembled_unmapped.fasta

    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Genome annotation runtime: $runtime min" >> $OUTPUT_DIR/${base}_time.log
fi


## Protein annotation -----------
if $PROT_ANN; then
    echo "====================================================================================="
    echo "Protein annotation with DIAMOND on unmapped sequences ..."
    start=`date +%s`
    mkdir -p dmnd_tmp
    diamond blastx --id 85 --query-cover 65 --min-score 60 \
        -p $n_threads -d $REF_DIR/nr \
        -q ${base}_contigs_unmapped.fasta \
        -o ${base}_contigs.dmdout \
        -f 6 -t dmnd_tmp -k 10


    diamond blastx --id 85 --query-cover 65 --min-score 60 \
        -p $n_threads -d $REF_DIR/nr \
        -q ${base}_unassembled_unmapped.fasta \
        -o ${base}_unassembled.dmdout \
        -f 6 -t dmnd_tmp -k 10 \

    $PYSCRIPT_DIR/7_Diamond_Protein_Map.py \
        $REF_DIR/nr \
        ${base}_contigs_map.tsv \
        ${base}_genes_map.tsv \
        ${base}_proteins.fasta \
        ${base}_contigs_unmapped.fasta \
        ${base}_contigs.dmdout \
        ${base}_contigs_unannotated.fasta \
        ${base}_unassembled_unmapped.fasta \
        ${base}_unassembled.dmdout \
        ${base}_unassembled_unannotated.fasta

    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Protein annotation (diamond) runtime: $runtime min" >> \
        $OUTPUT_DIR/${base}_time.log
fi

end0=`date +%s`
runtime0=$(((end0-start0)/60))
echo "======================================"
echo "ENTIRE WORKFLOW RUNTIME: $(runtime0) min" >> \
    $OUTPUT_DIR/${base}_time.log
