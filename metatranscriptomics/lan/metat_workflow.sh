#!/usr/bin/env bash
##
## Pipline for processing metatranscriptomics data based on workflow:
## https://github.com/ParkinsonLab/2017-Microbiome-Workshop
##
## author: nlhuong90@gmail.com
## date: 2/18/2018

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
export n_threads=${4:-4}
export use_sortmerna=true
export add_blat=false

## What Step to complete
export INDEX_DB=false
export TRIM=true
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
export PROT_ANN=false
export DIAMOND_REFSEQ=false
export DIAMOND_SEED=false

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
export SORTMERNA_DIR=$APP_DIR/sortmerna/build/Release/src/sortmerna/
export PYSCRIPT_DIR=$BASE_DIR/metatranscriptomics/pyscripts_edited

mkdir -p $DATA_DIR
mkdir -p $INPUT_DIR
mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR/QC

cd $OUTPUT_DIR

## Generate an index
if $INDEX_DB; then
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
fi


## Remove adapter seqs and trim low-quality seqs  ------------
if $TRIM && ! $paired; then
   ln -fs $APP_DIR/Trimmomatic-0.36/adapters/TruSeq3-SE.fa Adapters
   java -jar $APP_DIR/Trimmomatic-0.36/trimmomatic-0.36.jar \
       SE $INPUT_DIR/$input_fwd $OUTPUT_DIR/$fwd_trim.fastq \
       ILLUMINACLIP:Adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
   # Check reads quality
   $APP_DIR/FastQC/fastqc $INPUT_DIR/$input_fwd
   $APP_DIR/FastQC/fastqc $OUTPUT_DIR/$fwd_trim.fastq
   mv $OUTPUT_DIR/*.html $OUTPUT_DIR/QC/
   mv $OUTPUT_DIR/*.zip $OUTPUT_DIR/QC/
fi


## Do the same for paired ends ------------
if $TRIM && $paired; then
    ln -fs $APP_DIR/Trimmomatic-0.36/adapters/TruSeq3-SE.fa Adapters
    java -jar $APP_DIR/Trimmomatic-0.36/trimmomatic-0.36.jar \
        PE $INPUT_DIR/$input_fwd $INPUT_DIR/$input_rev \
        $OUTPUT_DIR/$fwd_paired_trim.fq $OUTPUT_DIR/$fw_unpaired_trim.fq \
        $OUTPUT_DIR/$rev_paired_trim.fq $OUTPUT_DIR/$rev_unpaired_trim.fq \
        ILLUMINACLIP:Adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
    # check reads quality
    $APP_DIR/FastQC/fastqc $INPUT_DIR/$input_fwd
    $APP_DIR/FastQC/fastqc $INPUT_DIR/$input_rev
    $APP_DIR/FastQC/fastqc $OUTPUT_DIR/$fwd_paired_trim.fq
    $APP_DIR/FastQC/fastqc $OUTPUT_DIR/$rev_paired_trim.fq
    mv $APP_DIR/*.html $OUTPUT_DIR/QC/
    mv $OUTPUT_DIR/*.zip $OUTPUT_DIR/QC/
fi


## Merge pairs ------------
if $MERGE_PAIRS; then
    $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
        --fastq_mergepairs $OUTPUT_DIR/$fwd_paired_trim.fq \
        --reverse $OUTPUT_DIR/$rev_paired_trim.fq \
        --fastqout $OUTPUT_DIR/$base_trim.fq \
        --fastqout_notmerged_fwd $OUTPUT_DIR/$base_unmerged_trim_fwd.fq \
        --fastqout_notmerged_rev $OUTPUT_DIR/$base_unmerged_trim_rev.fq \

    $APP_DIR/FastQC/fastqc $OUTPUT_DIR/$base_trim.fq
    mv $OUTPUT_DIR/*.html $OUTPUT_DIR/QC/
    mv $OUTPUT_DIR/*.zip $OUTPUT_DIR/QC/
fi


## Global quality filtering ------------
if $QUAL_FLTR; then
    $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
       --fastq_filter $OUTPUT_DIR/$base_trim.fq \
       --fastq_maxee 2.0 \
       --fastqout $OUTPUT_DIR/$base_qual.fq
fi


## Remove duplicate reads ------------
if $RM_DUPL; then
    $APP_DIR/cdhit/cd-hit-auxtools/cd-hit-dup \
        -i $OUTPUT_DIR/$base_qual.fq \
        -o $OUTPUT_DIR/$base_unique.fq
fi


## Remove unwanted vector contamination ------------
if $RM_VECTOR; then
    # align reads with vector db  and filter any reads that align to it
    bwa mem -t $n_threads $REF_DIR/UniVec_Core \
        $base_unique.fq > $base_univec_bwa.sam
    samtools view -bS $base_univec_bwa.sam > $base_univec_bwa.bam
    samtools fastq -n -F 4 -0 $base_univec_bwa_contam.fq $base_univec_bwa.bam
    samtools fastq -n -f 4 -0 $base_univec_bwa.fq $base_univec_bwa.bam

    # additional alignments for the reads with BLAT to filter out any remaining
    # reads that align to our vector contamination database but BLAT only accepts
    # fasta files so we have to convert our reads from fastq to fasta w/ vsearch
    $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
        --fastq_filter $base_univec_bwa.fq \
        --fastaout $base_univec_bwa.fasta

    blat $REF_DIR/UniVec_Core \
        $base_univec_bwa.fasta \
        -noHead -minIdentity=90 -minScore=65 \
        -fine -q=rna -t=dna -out=blast8 \
        $base_univec.blatout

    $PYSCRIPT_DIR/1_BLAT_Filter.py \
        $base_univec_bwa.fq \
        $base_univec.blatout \
        $base_univec_blat.fq \
        $base_univec_blat_contam.fq
fi

## Remove host reads ----------
# the same logic as the previous step
if $RM_HOST; then
    bwa mem -t $n_threads $REF_DIR/human_cds.fa \
        $base_univec_blat.fq > $base_human_bwa.sam
    samtools view -bS $base_human_bwa.sam > $base_human_bwa.bam
    samtools fastq -n -F 4 -0 $base_human_bwa_contam.fq $base_human_bwa.bam
    samtools fastq -n -f 4 -0 $base_human_bwa.fq $base_human_bwa.bam

    $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
        --fastq_filter $base_human_bwa.fq \
        --fastaout $base_human_bwa.fasta

    blat $REF_DIR/human_cds.fa \
        $base_human_bwa.fasta \
        -noHead -minIdentity=90 -minScore=65 \
        -fine -q=rna -t=dna -out=blast8 \
        $base_human.blatout

    $PYSCRIPT_DIR/1_BLAT_Filter.py \
        $base_human_bwa.fq \
        $base_human.blatout \
        $base_human_blat.fq \
        $base_human_blat_contam.fq
fi

## Remove rRNA seqs -----------
if $RM_rRNA %% $use_sortmerna; then
    $SORTMERNA_DIR/sortmerna \
    --ref $SORTMERNA_DIR/rRNA_databases/silva-bac-16s-id90.fasta:$SORTMERNA_DIR/index/silva-bac-16s-db \
        --reads $base_human_blat.fq \
        --aligned $base_unique_rRNA.fq \
        --other $base_unique_mRNA.fq \
        --fastx --num_alignments 0 --log -v
fi

if $RM_rRNA && ! $use_sortmerna; then
    $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
        --fastq_filter $base_human_blat.fq \
        --fastaout $base_human_blat.fasta

    $APP_DIR/infernal-1.1.2-linux-intel-gcc/binaries/cmsearch \
        -o $base_rRNA.log \
        --tblout $base_rRNA.infernalout \
        --anytrunc --rfam -E 0.001 \
        $REF_DIR/Rfam.cm \
        $base_human_blat.fasta

    $PYSCRIPT_DIR/2_Infernal_Filter.py \
        $base_human_blat.fq \
        $base_rRNA.infernalout \
        $base_unique_mRNA.fq \
        $base_unique_rRNA.fq
fi


## Rereplication (add back the repeated reads) ----------
if $REREPLICATION; then
    $PYSCRIPT_DIR/3_Reduplicate.py \
        $base_qual.fq \
        $base_unique_mRNA.fq \
        $base_unique.fastq.clstr \
        $base_mRNA.fq

    $APP_DIR/FastQC/fastqc $base_mRNA.fq
    mv $OUTPUT_DIR/*.html $OUTPUT_DIR/QC/
    mv $OUTPUT_DIR/*.zip $OUTPUT_DIR/QC/
fi

## Taxonomic Classification ------------
if $TAX_CLASS; then
   $APP_DIR/kaiju/bin/kaiju \
       -t $KAIJUBD_DIR/nodes.dmp \
       -f $KAIJUBD_DIR/kaiju_db.fmi \
       -i $base_mRNA.fq \
       -z $n_threads \
       -o $base_tax_class.tsv

   $PYSCRIPT_DIR/4_Constrain_Classification.py \
       genus \
       $base_tax_class.tsv \
       $KAIJUBD_DIR/nodes.dmp \
       $KAIJUBD_DIR/names.dmp \
       $base_genus_class.tsv

   $APP_DIR/kaiju/bin/kaijuReport \
       -t $KAIJUBD_DIR/nodes.dmp \
       -n $KAIJUBD_DIR/names.dmp \
       -i $base_genus_class.tsv \
       -o $base_genus_class_summary.txt \
       -r genus
fi

## Assemble reads into contigs -----------
if $ASSEMBLE; then
    # Build contigs
    $APP_DIR/SPAdes-3.11.1-Linux/bin/spades.py --rna \
        -s $base_mRNA.fq \
        -o $base_spades
    mv $base_spades/transcripts.fasta $base_contigs.fasta

    # Use the contigs as the database and align unassembled mRNA reads to it
    bwa index -a bwtsw $base_contigs.fasta
    bwa mem -t $n_threads $base_contigs.fasta $base_mRNA.fq > $base_contigs.sam

    # Make reads to contigs map
    $PYSCRIPT_DIR/5_Contig_Map.py \
        $base_mRNA.fq \
        $base_contigs.sam \
        $base_unassembled.fq \
        $base_contigs_map.tsv
fi


## Genome annotation -----------
if $GENOME_ANN; then
     # BWA utilizes nucleotide searches, we rely on a microbial genome database
    # Search DB
    bwa mem -t $n_threads $REF_DIR/microbial_all_cds.fasta \
        $base_contigs.fasta > $base_contigs_annotation_bwa.sam
    bwa mem -t $n_threads $REF_DIR/microbial_all_cds.fasta \
        $base_unassembled.fq > $base_unassembled_annotation_bwa.sam

    # IS THE FOLLOWING EVEN USEFUL? HOW TO ADD TO 6_BWA_Gene_Map.py ???
    if $add_blat; then
        # If we want to refine bwa alignments using blat
        # Extract reads that did not map to microbial genomes
        samtools view -bS \
           $base_contigs_annotation_bwa.sam > $base_contigs_annotation_bwa.bam
        samtools fastq -n -F 4 -0 \
            $base_contigs_annotation_aligned_bwa.fq \
            $base_contigs_annotation_bwa.bam
        samtools fastq -n -f 4 -0 \
            $base_contigs_annotation_bwa.fq \
            $base_contigs_annotation_bwa.bam

        $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
            --fastq_filter $base_contigs_annotation_bwa.fq \
            --fastaout $base_contigs_annotation_bwa.fasta

        # Additional alignments to reads not mapped by bwa
        blat $REF_DIR/microbial_all_cds.fasta \
            $base_contigs_annotation_bwa.fasta \
            -fine -q=rna -t=dna -out=blast8 \
            -noHead -minIdentity=90 -minScore=65  \
            $base_contigs_annotation_bwa.blatout

        $PYSCRIPT_DIR/1_BLAT_Filter.py \
            $base_contigs_annotation_bwa.fq \
            $base_contigs_annotation_bwa.blatout \
            $base_contigs_annotation_blat.fq \
            $base_contigs_annotation_aligned_blat.fq

        cat $base_contigs_annotation_aligned_bwa.fq \
            $base_contigs_annotation_aligned_blat.fq > \
            $base_contigs_annotation_aligned_bwa_blat.fq
        ## DO THE SAME FOR UNASSEMBLED???
    fi

    $PYSCRIPT_DIR/6_BWA_Gene_Map.py \
        $REF_DIR/microbial_all_cds.fasta \
        $base_contigs_map.tsv \
        $base_genes_map.tsv \
        $base_genes.fasta \
        $base_contigs.fasta \
        $base_contigs_annotation_bwa.sam \
        $base_contigs_unmapped.fasta \
        $base_unassembled.fq \
        $base_unassembled_annotation_bwa.sam \
        $base_unassembled_unmapped.fasta
fi


## Protein annotation -----------
if $PROT_ANN; then
    mkdir -p dmnd_tmp
    diamond blastx --id 85 --query-cover 65 --min-score 60 \
        -p $n_threads -d $REF_DIR/nr \
        -q $base_contigs_unmapped.fasta \
        -o $base_contigs.dmdout \
        -f 6 -t dmnd_tmp -k 10


    diamond blastx --id 85 --query-cover 65 --min-score 60 \
        -p $n_threads -d $REF_DIR/nr \
        -q $base_unassembled_unmapped.fasta \
        -o $base_unassembled.dmdout \
        -f 6 -t dmnd_tmp -k 10 \

    $PYSCRIPT_DIR/7_Diamond_Protein_Map.py \
        $REF_DIR/nr \
        $base_contigs_map.tsv \
        $base_genes_map.tsv \
        $base_proteins.fasta \
        $base_contigs_unmapped.fasta \
        $base_contigs.dmdout \
        $base_contigs_unannotated.fasta \
        $base_unassembled_unmapped.fasta \
        $base_unassembled.dmdout \
        $base_unassembled_unannotated.fasta
fi

## Additional annotation -----------
if $DIAMOND_REFSEQ; then
    # Align with RefSeq database
    mkdir -p dmnd_tmp
    diamond blastx -p $n_threads -db $REF_DIR/RefSeq_bac \
        -q $base_mRNA.fq -a $base.RefSeq \
        -t ./dmnd_tmp -k 1 --sensitive
    diamond view --daa $base.RefSeq.daa -f 6 \
        -o $base_refseq.dmdout
fi


if $DIAMOND_SEED; then
    # Align with SEED subsystem database
    mkdir -p dmnd_tmp
    diamond blastx -p $n_threads -db $REF_DIR/subsys_d \
        -q $base_mRNA.fq -a $base.Subsys \
        -t ./dmnd_tmp -k 1 --sensitive
    diamond view --daa $base.Subsys.daa -f 6 \
        -o $base_seed.dmdout
fi


