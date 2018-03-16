#!/usr/bin/env bash
##
## Pipline for processing metatranscriptomics data based on workflow:
## https://github.com/ParkinsonLab/2017-Microbiome-Workshop
##
## author: nlhuong90@gmail.com
## date: 2/18/2018

## DIRECTORIES ----------
if [ -z ${SCRATCH+x} ]; then
    #the variable $SCRATCH is unset
    echo Working on ICME cluster
    module load gcc/gcc6
    BASE_DIR=~/Projects/perturbation_16s
    APP_DIR=~/.local/bin/
    # PYSCRIPT_DIR=$BASE_DIR/metatranscriptomics/pyscripts_edited
else
    echo Working on SHERLOCK cluster
    BASE_DIR=$SCRATCH/Projects/perturbation_16s
    APP_DIR=$SCRATCH/applications/bin/
    # PYSCRIPT_DIR=$BASE_DIR/metatranscriptomics/pyscripts
    if [ $SHERLOCK == "1" ]; then
        CDHIT_DIR=$APP_DIR/for_shelock1/cdhit
        SORTMERNA_DIR=$APP_DIR/for_sherlock1/sortmerna
        KAIJU_DIR=$APP_DIR/for_sherlock1/kaiju/bin
        DIAMOND=$APP_DIR/for_sherlock1/diamond
    else 
        CDHIT_DIR=$APP_DIR/cdhit
        SORTMERNA_DIR=$APP_DIR/sortmerna
        KAIJU_DIR=$APP_DIR/kaiju/bin
        DIAMOND=$APP_DIR/diamond
    fi
fi

# Reference directories
REF_DIR=$BASE_DIR/data/databases
KAIJUBD_DB=$REF_DIR/kaijudb
# Python scripts directory
PYSCRIPT_DIR=$BASE_DIR/metatranscriptomics/pyscripts_edited
# SAMSA2 directory
SAMSA2=$APP_DIR/samsa2/python_scripts/
###############################################################################
## DEFAULTS (to be removed) ----------
INPUT_DIR=$BASE_DIR/data/metatranscriptomics/resilience/input/DBUr_Sub
OUTPUT_DIR=$BASE_DIR/data/metatranscriptomics/resilience/output/DBUr_Sub
# Fwd/Rev files
input_fwd=M3303_DBUsw_2r_TrM31_1P.fq.gz
input_rev=M3303_DBUsw_2r_TrM31_2P.fq.gz
###############################################################################

## PARAMETERS ----------
n_threads=5
rna_samples=true
use_sortmerna=true
extra_blat=false
index_db=false

## STEPS TO RUN ----------
TRIM=true
MERGE_PAIRS=true
QUAL_FLTR=true
RM_DUPL=true
RM_VECTOR=true
RM_HOST=true
RM_rRNA=true
REREPLICATION=true
TAX_CLASS=true
ASSEMBLE=true
GENOME_ANN=true
PROT_ANN=true
DIAMOND_REFSEQ=true
DIAMOND_SEED=true
DIAMOND_COUNT=true

## HELP DOCS ----------
usage() {
s=" "
tab="    "
echo
echo This program is a pipeline for processing raw metatranscriptomic reads.
echo
echo  "$s" -h --help  prints documentation
echo
echo  "$s" -i, --in_dir "$tab" location of input raw read files \(.fastq, .fq, .fasta, .fa\)
echo  "$s" -o, --out_dir "$tab" directory of output files
echo  "$s" -f, --fwd "$tab" forward reads file
echo  "$s" -r, --rev "$tab" \(optional\) reverse reads file
echo
echo Additional options:
echo
echo  "$s" -t X "$tab" set number of parallel threads for the pipeline X \(default:5\)
echo  "$s" --metagenome "$tab" whether to input files are metagenomic, NOT metatranscriptomic reads \(default:false\).
echo  "$s" --no-assembly "$tab" whether to apply DIAMOND directly to processed short reads w/o assembling first \(default:false\).
echo  "$s" --index-db "$tab" whether to index databases \(default:false should be done ahead\).
echo  "$s" --sortmerna "$tab" use SortMeRNA instead of Infernal for ribodepletion step \(default:true\).
echo  "$s" --extra-blat "$tab" use extra BLAT step during genome annotation \(default:false\).
echo
}

## ARGUMENTS ----------
#POSITIONAL=()
while [[ "$#" -gt 0 ]]
do
key="$1"
case $key in
    -h|-\?|--help)
    usage
    exit 1
    ;;
    -i|--in_dir)
    INPUT_DIR="$2"
    if [ ! -d "$INPUT_DIR" ]; then
       echo "Input directory $INPUT_DIR does not exist."
       exit 1
    fi
    shift # past argument
    shift # past value
    ;;
    -o|--out_dir)
    OUTPUT_DIR="$2"
    shift # past argument
    shift # past value
    ;;
    -f|--fwd)
    input_fwd="$2"
    if [ ! -f "$INPUT_DIR/$input_fwd" ]; then
       echo "Forward file $INPUT_DIR/$input_fwd does not exist."
       exit 1
    fi
    shift # past argument
    shift # past value
    ;;
    -r|--rev)
    input_rev="$2"
    if [ ! -f $INPUT_DIR/$input_rev ]; then
       echo "Reverse file $INPUT_DIR/$input_rev does not exist."
       exit 1
    fi
    shift # past argument
    shift # past value
    ;;
    -t|--nThreads)
    n_threads="$2"
    if [ ! -n "$n_threads" ]; then
        printf 'ERROR: "-t" requires a non-empty integer argument.\n' >&2
        usage
        exit 1
    fi
    shift # past argument
    shift # past value
    ;;
    --metagenome)
    rna_samples=false
    RM_rRNA=false
    PROT_ANN=false
    shift # past argument
    ;;
    --no-assembly)
    ASSEMBLE=false
    GENOME_ANN=false
    PROT_ANN=false
    DIAMOND_REFSEQ=true
    DIAMOND_SEED=true
    shift # past argument
    ;;
    --sortmerna)
    use_sortmerna=true
    shift # past argument
    ;;
    --index-db)
    index_db=true
    shift # past argument
    ;;
    --extra-blat)
    extra_blat=true
    shift # past argument
    ;;
    -*) echo "unknown option: $1" >&2; exit 1;;
    *) handle_argument "$1"; shift 1;;
    #POSITIONAL+=("$1") # save it in an array for later
    #shift # past argument
    #;;
esac
done
#set -- "${POSITIONAL[@]}" # restore positional parameters

if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$input_fwd" ] ||  [ -z "$input_rev" ];
then
    echo "Missing arguments: workflow.sh requires -i, -o, -f, -r arguments specification."
fi

###############################################################################

base=${input_fwd%_1P.fq.gz}
fwd=${input_fwd%.fq.gz}
rev=${input_rev%.fq.gz}

if [ $input_rev = "" ]; then
    paired=false
    base=fwd
else
    paired=true
fi


echo STARTING DATA PROCESSING PIPELINE
echo " "
echo =======================================================================
echo "   "
echo PAIRED READS     = "${paired}"
echo FWD FILE         = "${input_fwd}"
echo REV FILE         = "${input_rev}"
echo NO. THREADS      = "${n_threads}"
echo " "

## Make output direcories ----------------
mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR/time
mkdir -p $OUTPUT_DIR/main
mkdir -p $OUTPUT_DIR/QC
mkdir -p $OUTPUT_DIR/trimmed/
mkdir -p $OUTPUT_DIR/aligned/
mkdir -p $OUTPUT_DIR/dmnd_tmp/
mkdir -p $OUTPUT_DIR/counts/

if [ $GENOME_ANN ]; then mkdir -p $OUTPUT_DIR/genome/ ; fi
if [ $PROT_ANN ] || [ $DIAMOND_REFSEQ ] || [ $DIAMOND_SEED ]; then
     mkdir -p $OUTPUT_DIR/diamond/
fi
if [ $PROT_ANN ]; then mkdir -p $OUTPUT_DIR/protein ; fi

cd $OUTPUT_DIR

###############################################################################
if [ $rna_samples ]; then
    echo "RUNTIME FOR METATRANSCTIPTOMIC PIPELINE:" >> \
        $OUTPUT_DIR/time/${base}_time.log
else
    echo "RUNTIME FOR METAGENOMIC PIPELINE:" >> \
        $OUTPUT_DIR/time/${base}_time.log
fi
echo "FOR sample: $base" >> $OUTPUT_DIR/time/${base}_time.log
echo ========================================== >> $OUTPUT_DIR/time/${base}_time.log


start0=`date +%s`
## Generate an index
if $index_db; then
    echo =======================================================================
    echo Indexing databases ...
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
    makeblastdb -in $REF_DIR/microbial_all_cds.fasta -dbtype nucl
    # Build diamond compatible rRNA_databases
    diamond makedb -p $n_threads \
        --in $REF_DIR/RefSeq_bac.fa \
        --db $REF_DIR/RefSeq_bac
    diamond makedb -p $n_threads \
        --in $REF_DIR/subsys_db.fa \
        --db $REF_DIR/subsys_db
    diamond makedb -p $n_threads \
        --in $REF_DIR/nr --db $REF_DIR/nr
    # diamond makedb -p $n_threads \
    #     --in $REF_DIR/uniref100.fasta \
    #     --db $REF_DIR/uniref100
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Reference database indexing: $runtime min" >> \
        $OUTPUT_DIR/time/${base}_time.log
fi


## Remove adapter seqs and trim low-quality seqs  ------------
if $TRIM && ! $paired && [ ! -f $OUTPUT_DIR/trimmed/${fwd}_trim.fastq ]; then
   echo =======================================================================
   echo Trimming and removing adapters ...
   start=`date +%s`
   ln -fs $APP_DIR/Trimmomatic-0.36/adapters/TruSeq3-SE.fa Adapters
   java -jar $APP_DIR/Trimmomatic-0.36/trimmomatic-0.36.jar \
       SE $INPUT_DIR/$input_fwd $OUTPUT_DIR/trimmed/${fwd}_trim.fastq \
       ILLUMINACLIP:Adapters:2:30:10 LEADING:3 \
       TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
   # Check reads quality
   $APP_DIR/FastQC/fastqc $INPUT_DIR/$input_fwd
   $APP_DIR/FastQC/fastqc $OUTPUT_DIR/trimmed/${fwd}_trim.fastq
   mv $INPUT_DIR/*.html $OUTPUT_DIR/QC/
   mv $INPUT_DIR/*.zip $OUTPUT_DIR/QC/
   mv $OUTPUT_DIR/trimmed/*.html $OUTPUT_DIR/QC/
   mv $OUTPUT_DIR/trimmed/*.zip $OUTPUT_DIR/QC/
   end=`date +%s`
   runtime=$(((end-start)/60))
   echo "Trimming and QC: $runtime min" >> $OUTPUT_DIR/time/${base}_time.log
   echo ========================================== >> \
       $OUTPUT_DIR/time/${base}_time.log
fi


## Do the same for paired ends ------------
if $TRIM && $paired && [ ! -f $OUTPUT_DIR/trimmed/${rev}_paired_trim.fq ]; then
    echo =======================================================================
    echo Trimming and removing adapters ...
    start=`date +%s`
    ln -fs $APP_DIR/Trimmomatic-0.36/adapters/TruSeq3-SE.fa Adapters
    java -jar $APP_DIR/Trimmomatic-0.36/trimmomatic-0.36.jar \
        PE $INPUT_DIR/$input_fwd $INPUT_DIR/$input_rev \
        $OUTPUT_DIR/trimmed/${fwd}_paired_trim.fq \
        $OUTPUT_DIR/trimmed/${fwd}_unpaired_trim.fq \
        $OUTPUT_DIR/trimmed/${rev}_paired_trim.fq \
        $OUTPUT_DIR/trimmed/${rev}_unpaired_trim.fq \
        ILLUMINACLIP:Adapters:2:30:10 LEADING:3 \
        TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
    # check reads quality
    $APP_DIR/FastQC/fastqc $INPUT_DIR/$input_fwd
    $APP_DIR/FastQC/fastqc $INPUT_DIR/$input_rev
    $APP_DIR/FastQC/fastqc $OUTPUT_DIR/trimmed/${fwd}_paired_trim.fq
    $APP_DIR/FastQC/fastqc $OUTPUT_DIR/trimmed/${rev}_paired_trim.fq
    mv $INPUT_DIR/*.html $OUTPUT_DIR/QC/
    mv $INPUT_DIR/*.zip $OUTPUT_DIR/QC/
    mv $OUTPUT_DIR/trimmed/*.html $OUTPUT_DIR/QC/
    mv $OUTPUT_DIR/trimmed/*.zip $OUTPUT_DIR/QC/
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Trimming and QC: $runtime min" >> $OUTPUT_DIR/time/${base}_time.log
    echo ========================================== >> \
        $OUTPUT_DIR/time/${base}_time.log
fi


## Merge pairs ------------
if $MERGE_PAIRS && [ ! -f $OUTPUT_DIR/trimmed/${base}_unmerged_trim_rev.fq ]; then
    echo =======================================================================
    echo Merging paired reads ...
    start=`date +%s`
    $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
        --fastq_mergepairs $OUTPUT_DIR/trimmed/${fwd}_paired_trim.fq \
        --reverse $OUTPUT_DIR/trimmed/${rev}_paired_trim.fq \
        --fastqout $OUTPUT_DIR/trimmed/${base}_trim.fq \
        --fastqout_notmerged_fwd $OUTPUT_DIR/trimmed/${base}_unmerged_trim_fwd.fq \
        --fastqout_notmerged_rev $OUTPUT_DIR/trimmed/${base}_unmerged_trim_rev.fq \

    $APP_DIR/FastQC/fastqc $OUTPUT_DIR/trimmed/${base}_trim.fq
    mv $OUTPUT_DIR/trimmed/*.html $OUTPUT_DIR/QC/
    mv $OUTPUT_DIR/trimmed/*.zip $OUTPUT_DIR/QC/
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Merging reads: $runtime min" >> $OUTPUT_DIR/time/${base}_time.log
    echo ========================================== >> \
        $OUTPUT_DIR/time/${base}_time.log
fi


## Global quality filtering ------------
if $QUAL_FLTR && [ ! -f $OUTPUT_DIR/main/${base}_qual.fq ]; then
    echo =======================================================================
    echo Quality filtering ...
    start=`date +%s`
    $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
       --fastq_filter $OUTPUT_DIR/trimmed/${base}_trim.fq \
       --fastq_maxee 2.0 \
       --fastqout $OUTPUT_DIR/main/${base}_qual.fq
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Quality filtering: $runtime min" >> \
        $OUTPUT_DIR/time/${base}_time.log
    echo ========================================== >> \
        $OUTPUT_DIR/time/${base}_time.log
fi


## Remove duplicate reads ------------
if $RM_DUPL && [ ! -f $OUTPUT_DIR/main/${base}_unique.fq ]; then
    echo =======================================================================
    echo Remove duplicates ...
    start=`date +%s`
    $CDHIT_DIR/cd-hit-auxtools/cd-hit-dup \
        -i $OUTPUT_DIR/main/${base}_qual.fq \
        -o $OUTPUT_DIR/main/${base}_unique.fq
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Remove duplicates: $runtime min" >> \
        $OUTPUT_DIR/time/${base}_time.log
    echo ========================================== >> \
        $OUTPUT_DIR/time/${base}_time.log
fi


## Remove unwanted vector contamination ------------
if $RM_VECTOR && [ ! -f $OUTPUT_DIR/aligned/${base}_univec_blat.fq ]; then
    echo =======================================================================
    echo Remove vector sequences ...
    start=`date +%s`
    # align reads with vector db  and filter any reads that align to it
    bwa mem -t $n_threads $REF_DIR/UniVec_Core \
        $OUTPUT_DIR/main/${base}_unique.fq > \
        $OUTPUT_DIR/aligned/${base}_univec_bwa.sam
    samtools view -bS \
        $OUTPUT_DIR/aligned/${base}_univec_bwa.sam > \
        $OUTPUT_DIR/aligned/${base}_univec_bwa.bam
    samtools fastq -n -F 4 -0 \
        $OUTPUT_DIR/aligned/${base}_univec_bwa_contam.fq \
        $OUTPUT_DIR/aligned/${base}_univec_bwa.bam
    samtools fastq -n -f 4 -0 \
        $OUTPUT_DIR/aligned/${base}_univec_bwa.fq \
        $OUTPUT_DIR/aligned/${base}_univec_bwa.bam

    # additional alignments for the reads with BLAT to filter out any remaining
    # reads that align to our vector contamination database but BLAT only
    # accepts .fasta so we convert reads from fastq to fasta w/ vsearch
    $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
        --fastq_filter $OUTPUT_DIR/aligned/${base}_univec_bwa.fq \
        --fastaout $OUTPUT_DIR/aligned/${base}_univec_bwa.fasta

    blat $REF_DIR/UniVec_Core \
        $OUTPUT_DIR/aligned/${base}_univec_bwa.fasta \
        -noHead -minIdentity=90 -minScore=65 \
        -fine -q=rna -t=dna -out=blast8 \
        $OUTPUT_DIR/aligned/${base}_univec.blatout

    $PYSCRIPT_DIR/1_BLAT_Filter.py \
        $OUTPUT_DIR/aligned/${base}_univec_bwa.fq \
        $OUTPUT_DIR/aligned/${base}_univec.blatout \
        $OUTPUT_DIR/aligned/${base}_univec_blat.fq \
        $OUTPUT_DIR/aligned/${base}_univec_blat_contam.fq

    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Removing vectors: $runtime min" >> \
        $OUTPUT_DIR/time/${base}_time.log
    echo ========================================== >> \
        $OUTPUT_DIR/time/${base}_time.log
fi

## Remove host reads ----------
# the same logic as the previous step
if $RM_HOST && [ ! -f $OUTPUT_DIR/aligned/${base}_human_blat.fq ]; then
    echo =======================================================================
    echo Remove host sequences ...
    start=`date +%s`
    bwa mem -t $n_threads $REF_DIR/human_cds.fa \
        $OUTPUT_DIR/aligned/${base}_univec_blat.fq > \
        $OUTPUT_DIR/aligned/${base}_human_bwa.sam
    samtools view -bS \
        $OUTPUT_DIR/aligned/${base}_human_bwa.sam > \
        $OUTPUT_DIR/aligned/${base}_human_bwa.bam
    samtools fastq -n -F 4 -0 \
        $OUTPUT_DIR/aligned/${base}_human_bwa_contam.fq \
        $OUTPUT_DIR/aligned/${base}_human_bwa.bam
    samtools fastq -n -f 4 -0 \
        $OUTPUT_DIR/aligned/${base}_human_bwa.fq \
        $OUTPUT_DIR/aligned/${base}_human_bwa.bam

    $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
        --fastq_filter $OUTPUT_DIR/aligned/${base}_human_bwa.fq \
        --fastaout $OUTPUT_DIR/aligned/${base}_human_bwa.fasta

    blat $REF_DIR/human_cds.fa \
        $OUTPUT_DIR/aligned/${base}_human_bwa.fasta \
        -noHead -minIdentity=90 -minScore=65 \
        -fine -q=rna -t=dna -out=blast8 \
        $OUTPUT_DIR/aligned/${base}_human.blatout

    $PYSCRIPT_DIR/1_BLAT_Filter.py \
        $OUTPUT_DIR/aligned/${base}_human_bwa.fq \
        $OUTPUT_DIR/aligned/${base}_human.blatout \
        $OUTPUT_DIR/aligned/${base}_human_blat.fq \
        $OUTPUT_DIR/aligned/${base}_human_blat_contam.fq

    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Removing host: $runtime min" >> $OUTPUT_DIR/time/${base}_time.log
    echo ========================================== >> \
        $OUTPUT_DIR/time/${base}_time.log
fi

## Remove rRNA seqs -----------
if $RM_rRNA && $use_sortmerna && [ ! -f $OUTPUT_DIR/main/${base}_unique_mRNA.fq ]; then
    echo =======================================================================
    echo Ribodepletion with SortMeRNA ...
    start=`date +%s`
    $SORTMERNA_DIR/build/Release/src/sortmerna/sortmerna \
        --ref $SORTMERNA_DIR/rRNA_databases/silva-bac-16s-id90.fasta,$SORTMERNA_DIR/index/silva-bac-16s \
        --reads $OUTPUT_DIR/aligned/${base}_human_blat.fq \
        --aligned $OUTPUT_DIR/main/${base}_unique_rRNA \
        --other $OUTPUT_DIR/main/${base}_unique_mRNA \
        --fastx --num_alignments 0 --log \
        -a $n_threads -v
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "SortMeRNA ribodepletion: $runtime min" >> \
        $OUTPUT_DIR/time/${base}_time.log
    echo ========================================== >> \
        $OUTPUT_DIR/time/${base}_time.log
fi

if $RM_rRNA && ! $use_sortmerna && [ ! -f $OUTPUT_DIR/main/${base}_unique_mRNA.fq ]; then
    echo =======================================================================
    echo Ribodepletion with infernalout ...
    start=`date +%s`
    mkdir -p $OUTPUT_DIR/infernal
    $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
        --fastq_filter $OUTPUT_DIR/aligned/${base}_human_blat.fq \
        --fastaout $OUTPUT_DIR/aligned/${base}_human_blat.fasta

    $APP_DIR/infernal-1.1.2-linux-intel-gcc/binaries/cmsearch \
        -o $OUTPUT_DIR/infernal/${base}_rRNA.log \
        --tblout $OUTPUT_DIR/infernal/${base}_rRNA.infernalout \
        --anytrunc --rfam -E 0.001 --cpu $n_threads \
        $REF_DIR/Rfam.cm \
        $OUTPUT_DIR/aligned/${base}_human_blat.fasta

    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Infernal ribodepletion: $runtime min" >> \
        $OUTPUT_DIR/time/${base}_time.log

    start=`date +%s`
    $PYSCRIPT_DIR/2_Infernal_Filter.py \
        $OUTPUT_DIR/aligned/${base}_human_blat.fq \
        $OUTPUT_DIR/infernal/${base}_rRNA.infernalout \
        $OUTPUT_DIR/main/${base}_unique_mRNA.fq \
        $OUTPUT_DIR/main/${base}_unique_rRNA.fq

    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Python infernal filter: $runtime min" >> \
        $OUTPUT_DIR/time/${base}_time.log
    echo ========================================== >> \
        $OUTPUT_DIR/time/${base}_time.log
fi


## Rereplication (add back the repeated reads) ----------
if $REREPLICATION; then
    if [ $rna_samples ]; then
        outfile=$OUTPUT_DIR/main/${base}_mRNA.fq
    else
        outfile=$OUTPUT_DIR/main/${base}_fltr.fq
    fi
    
    if [ ! -f $outfile ]; then
        echo =======================================================================
        echo Rereplication ...
        start=`date +%s`
        $PYSCRIPT_DIR/3_Reduplicate.py \
            $OUTPUT_DIR/main/${base}_qual.fq \
            $OUTPUT_DIR/aligned/${base}_human_blat.fq \
            $OUTPUT_DIR/main/${base}_unique.fq.clstr \
            $outfile

        $APP_DIR/FastQC/fastqc $outfile
        mv $OUTPUT_DIR/*.html $OUTPUT_DIR/QC/
        mv $OUTPUT_DIR/*.zip $OUTPUT_DIR/QC/
        end=`date +%s`
        runtime=$(((end-start)/60))
        echo "Rereplication: $runtime min" >> $OUTPUT_DIR/time/${base}_time.log
        echo ========================================== >> \
            $OUTPUT_DIR/time/${base}_time.log
    fi
fi


## Taxonomic Classification ------------
if $TAX_CLASS && [ ! -f $OUTPUT_DIR/taxonomy/${base}_genus_class_summary.txt ]; then
    mkdir -p taxonomy
    echo =======================================================================
    echo Taxonomy classification ...

    if [ $rna_samples ]; then
        infile=$OUTPUT_DIR/main/${base}_mRNA.fq
    else
        infile=$OUTPUT_DIR/main/${base}_fltr.fq
    fi

    start=`date +%s`
    $KAIJU_DIR/kaiju \
        -t $KAIJUBD_DB/nodes.dmp \
        -f $KAIJUBD_DB/kaiju_db.fmi \
        -i $infile \
        -z $n_threads \
        -o $OUTPUT_DIR/taxonomy/${base}_tax_class.tsv
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Kaiju genome annotation: $runtime min" >> $OUTPUT_DIR/time/${base}_time.log

    start=`date +%s`
    $PYSCRIPT_DIR/4_Constrain_Classification.py \
        genus \
        $OUTPUT_DIR/taxonomy/${base}_tax_class.tsv \
        $KAIJUBD_DB/nodes.dmp \
        $KAIJUBD_DB/names.dmp \
        $OUTPUT_DIR/taxonomy/${base}_genus_class.tsv
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Python script process taxonomy result: $runtime min" >> \
        $OUTPUT_DIR/time/${base}_time.log

    start=`date +%s`
    $KAIJU_DIR/kaijuReport \
        -t $KAIJUBD_DB/nodes.dmp \
        -n $KAIJUBD_DB/names.dmp \
        -i $OUTPUT_DIR/taxonomy/${base}_genus_class.tsv \
        -o $OUTPUT_DIR/taxonomy/${base}_genus_class_summary.txt \
        -r genus
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Aggregate taxonomy res: $runtime min" >> \
        $OUTPUT_DIR/time/${base}_time.log
    echo ========================================== >> \
        $OUTPUT_DIR/time/${base}_time.log
fi


## Assemble reads into contigs -----------
if $ASSEMBLE && [ ! -f $OUTPUT_DIR/assembled/${base}_contigs_map.tsv ]; then
    echo =======================================================================
    echo Contigs assembly ...

    # Build contigs
    start=`date +%s`
    if [ $rna_samples ]; then
        infile=$OUTPUT_DIR/main/${base}_mRNA.fq
        $APP_DIR/SPAdes-3.11.1-Linux/bin/spades.py \
            --rna -t $n_threads \
            -s $infile \
            -o $OUTPUT_DIR/assembled/${base}_spades
    else
        infile=$OUTPUT_DIR/assembled/${base}_spades
        $APP_DIR/SPAdes-3.11.1-Linux/bin/spades.py \
            -t $n_threads \
            -s $infile \
            -o $OUTPUT_DIR/assembled/${base}_spades
    fi
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "SPADES assembly: $runtime min" >> $OUTPUT_DIR/time/${base}_time.log

    mv $OUTPUT_DIR/assembled/${base}_spades/transcripts.fasta \
        $OUTPUT_DIR/assembled/${base}_contigs.fasta

    # Use the contigs as the database and align unassembled reads to it
    start=`date +%s`
    bwa index -a bwtsw $OUTPUT_DIR/assembled/${base}_contigs.fasta
    bwa mem -t $n_threads \
        $OUTPUT_DIR/assembled/${base}_contigs.fasta \
        $infile > $OUTPUT_DIR/assembled/${base}_contigs.sam
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Indexing contig and alingning reads: $runtime min" >> \
        $OUTPUT_DIR/time/${base}_time.log
    
    #start=`date +%s`
    #samtools view -bS \
    #    $OUTPUT_DIR/assembled/${base}_contigs.sam > \
    #    $OUTPUT_DIR/assembled/${base}_contigs.bam

    #samtools fastq -n -f 4 -0 \
    #    $OUTPUT_DIR/assembled/${base}_unassembled.fq \
    #    $OUTPUT_DIR/assembled/${base}_contigs.bam
    #end=`date +%s`
    #runtime=$(((end-start)/60))
    #echo "Samtools extract unassembled reads: $runtime min" >> \
    #    $OUTPUT_DIR/time/${base}_time.log
    #echo ========================================== >> \
    #    $OUTPUT_DIR/time/${base}_time.log
    
    # Make reads to contigs map
    start=`date +%s`
    $PYSCRIPT_DIR/5_Contig_Map.py \
        $infile \
        $OUTPUT_DIR/assembled/${base}_contigs.sam \
        $OUTPUT_DIR/assembled/${base}_unassembled.fq \
        $OUTPUT_DIR/assembled/${base}_contigs_map.tsv
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "Python aggregate assembly result: $runtime min" >> \
        $OUTPUT_DIR/time/${base}_time.log
    echo ========================================== >> \
        $OUTPUT_DIR/time/${base}_time.log
fi


## Genome annotation -----------
if $GENOME_ANN && [ ! -f $OUTPUT_DIR/assembled/${base}_unassembled_unmapped.fq ]; then
    echo =======================================================================
    echo Genome annotation with BWA ...
    mkdir -p $OUTPUT_DIR/genome/
    start=`date +%s`
    # BWA utilizes nucleotide searches, we rely on a microbial genome database
    # Search DB
    bwa mem -t $n_threads $REF_DIR/microbial_all_cds.fasta \
        $OUTPUT_DIR/assembled/${base}_contigs.fasta > \
        $OUTPUT_DIR/assembled/${base}_contigs_annotation_bwa.sam
    bwa mem -t $n_threads $REF_DIR/microbial_all_cds.fasta \
        $OUTPUT_DIR/assembled/${base}_unassembled.fq > \
        $OUTPUT_DIR/assembled/${base}_unassembled_annotation_bwa.sam
    samtools view -bS \
        $OUTPUT_DIR/assembled/${base}_contigs_annotation_bwa.sam > \
        $OUTPUT_DIR/assembled/${base}_contigs_annotation_bwa.bam
    samtools fastq -n -F 4 -0 \
        $OUTPUT_DIR/assembled/${base}_contigs_bwa_aligned.fq \
        $OUTPUT_DIR/assembled/${base}_contigs_annotation_bwa.bam
    samtools fastq -n -f 4 -0 \
        $OUTPUT_DIR/assembled/${base}_contigs_unmapped.fq \
        $OUTPUT_DIR/assembled/${base}_contigs_annotation_bwa.bam

    samtools view -bS \
        $OUTPUT_DIR/assembled/${base}_unassembled_annotation_bwa.sam > \
        $OUTPUT_DIR/assembled/${base}_unassembled_annotation_bwa.bam
    samtools fastq -n -F 4 -0 \
        $OUTPUT_DIR/assembled/${base}_unassembled_bwa_aligned.fq \
        $OUTPUT_DIR/assembled/${base}_unassembled_annotation_bwa.bam
    samtools fastq -n -f 4 -0 \
        $OUTPUT_DIR/assembled/${base}_unassembled_unmapped.fq \
        $OUTPUT_DIR/assembled/${base}_unassembled_annotation_bwa.bam

    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "BWA genome alignment: $runtime min" >> \
        $OUTPUT_DIR/time/${base}_time.log

    # start=`date +%s`
    # $PYSCRIPT_DIR/6_BWA_Gene_Map.py \
    #     $REF_DIR/microbial_all_cds.fasta \
    #     $OUTPUT_DIR/assembled/${base}_contigs_map.tsv \
    #     $OUTPUT_DIR/genome/${base}_genes_map.tsv \
    #     $OUTPUT_DIR/genome/${base}_genes.fasta \
    #     $OUTPUT_DIR/assembled/${base}_contigs.fasta \
    #     $OUTPUT_DIR/assembled/${base}_contigs_annotation_bwa.sam \
    #     $OUTPUT_DIR/assembled/${base}_contigs_unmapped.fasta \
    #     $OUTPUT_DIR/assembled/${base}_unassembled.fq \
    #     $OUTPUT_DIR/assembled/${base}_unassembled_annotation_bwa.sam \
    #     $OUTPUT_DIR/assembled/${base}_unassembled_unmapped.fasta
    # end=`date +%s`
    # runtime=$(((end-start)/60))
    # echo "Python aggregate genome alignment results: $runtime min" >> \
    #     $OUTPUT_DIR/time/${base}_time.log

    # IS THE FOLLOWING EVEN USEFUL? HOW TO ADD TO 6_BWA_Gene_Map.py ???
    if $extra_blat; then
        # If we want to refine bwa alignments using blat
        # Extract reads that did not map to microbial genomes
        start=`date +%s`
        blat $REF_DIR/microbial_all_cds.fasta \
            $OUTPUT_DIR/assembled/${base}_contigs_unmapped.fasta \
            -fine -q=rna -t=dna -out=blast8 \
            -noHead -minIdentity=90 -minScore=65  \
            $OUTPUT_DIR/assembled/${base}_contigs_unmapped_bwa.blatout

        blat $REF_DIR/microbial_all_cds.fasta \
            $OUTPUT_DIR/assembled/${base}_unassembled_unmapped.fasta \
            -fine -q=rna -t=dna -out=blast8 \
            -noHead -minIdentity=90 -minScore=65  \
            $OUTPUT_DIR/assembled/${base}_unassembled_unmapped_bwa.blatout

        $PYSCRIPT_DIR/1_BLAT_Filter.py \
            $OUTPUT_DIR/assembled/${base}_contigs_unmapped.fq \
            $OUTPUT_DIR/assembled/${base}_contigs_unmapped_bwa.blatout \
            $OUTPUT_DIR/assembled/${base}_contigs_unmapped_bwa_blat.fq \
            $OUTPUT_DIR/assembled/${base}_contigs_unmapped_bwa_mapped_blat.fq

        $PYSCRIPT_DIR/1_BLAT_Filter.py \
            $OUTPUT_DIR/assembled/${base}_unassembled_unmapped.fq \
            $OUTPUT_DIR/assembled/${base}_unassembled_unmapped_bwa.blatout \
            $OUTPUT_DIR/assembled/${base}_unassembled_unmapped_bwa_blat.fq \
            $OUTPUT_DIR/assembled/${base}_unassembled_unmapped_bwa_mapped_blat.fq

        rm $OUTPUT_DIR/assembled/${base}_contigs_unmapped.fasta
        rm $OUTPUT_DIR/assembled/${base}_unassembled_unmapped.fasta

        $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
            --fastq_filter $OUTPUT_DIR/assembled/${base}_contigs_unmapped_bwa_blat.fq \
            --fastaout $OUTPUT_DIR/assembled/${base}_contigs_unmapped.fasta

        $APP_DIR/vsearch-2.7.0-linux-x86_64/bin/vsearch \
            --fastq_filter $OUTPUT_DIR/assembled/${base}_unassembled_unmapped_bwa_blat.fq \
            --fastaout $OUTPUT_DIR/assembled/${base}_unassembled_unmapped.fasta

        end=`date +%s`
        runtime=$(((end-start)/60))
        echo "Extra BLAT genome alignment: $runtime min" >> \
            $OUTPUT_DIR/time/${base}_time.log
    fi
    echo ========================================== >> \
        $OUTPUT_DIR/time/${base}_time.log
fi


## Protein annotation -----------
if $PROT_ANN && [ ! -f $OUTPUT_DIR/diamond/${base}_nr_unassembled.dmdout ]; then
    echo =======================================================================
    echo Protein annotation with DIAMOND on unmapped sequences ...

    start=`date +%s`
    $DIAMOND blastx --id 85 --query-cover 65 --min-score 60 \
        --threads $n_threads \
        -d $REF_DIR/nr \
        -q $OUTPUT_DIR/assembled/${base}_contigs_unmapped.fq \
        -o $OUTPUT_DIR/diamond/${base}_nr_contigs.dmdout \
        -f 6 -t $OUTPUT_DIR/dmnd_tmp -k 10

    $DIAMOND blastx --id 85 --query-cover 65 --min-score 60 \
        --threads $n_threads \
        -d $REF_DIR/nr \
        -q $OUTPUT_DIR/assembled/${base}_unassembled_unmapped.fq \
        -o $OUTPUT_DIR/diamond/${base}_nr_unassembled.dmdout \
        -f 6 -t $OUTPUT_DIR/dmnd_tmp -k 10
    
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "DIAMOND (NR) protein alignment: $runtime min" >> \
        $OUTPUT_DIR/time/${base}_time.log

    # start=`date +%s`
    # $PYSCRIPT_DIR/7_Diamond_Protein_Map.py \
    #     $REF_DIR/nr \
    #     $OUTPUT_DIR/assembled/${base}_contigs_map.tsv \
    #     $OUTPUT_DIR/genome/${ibase}_genes_map.tsv \
    #     $OUTPUT_DIR/proteins/${base}_proteins.fasta \
    #     $OUTPUT_DIR/assembled/${base}_contigs_unmapped.fasta \
    #     $OUTPUT_DIR/diamond/${base}_contigs.dmdout \
    #     $OUTPUT_DIR/proteins/${base}_contigs_unannotated.fasta \
    #     $OUTPUT_DIR/assembled/${base}_unassembled_unmapped.fasta \
    #     $OUTPUT_DIR/diamond/${base}_unassembled.dmdout \
    #     $OUTPUT_DIR/proteins/${base}_unassembled_unannotated.fasta
    # end=`date +%s`
    # runtime=$(((end-start)/60))
    # echo "Python aggregating protein annotation results: $runtime min" >> \
    #     $OUTPUT_DIR/time/${base}_time.log

    echo ========================================== >> \
        $OUTPUT_DIR/time/${base}_time.log
fi


## DIAMOND annotation -----------
if $DIAMOND_REFSEQ && [ ! -f $OUTPUT_DIR/diamond/${base}_refseq.dmdout ]; then
    echo =======================================================================
    echo Refseq annotation with DIAMOND ...

    if [ $rna_samples ]; then
        infile=$OUTPUT_DIR/main/${base}_mRNA.fq
    else
        infile=$OUTPUT_DIR/main/${base}_fltr.fq
    fi

    start=`date +%s`
    # Align with RefSeq database
    $DIAMOND blastx --threads $n_threads \
        -d $REF_DIR/RefSeq_bac \
        -q $infile \
        -o $OUTPUT_DIR/diamond/${base}_refseq.dmdout \
        -t $OUTPUT_DIR/dmnd_tmp -f 6 -k 1 --sensitive

    # $DIAMOND blastx -p $n_threads -d $REF_DIR/RefSeq_bac \
    #     -q $infile \
    #     -a $OUTPUT_DIR/diamond/${base}.RefSeq \
    #     -t ./dmnd_tmp -k 1 --sensitive
    # $DIAMOND view --daa ${base}.RefSeq.daa -f 6 \
    #     -o $OUTPUT_DIR/diamond/${base}_refseq.dmdout

    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "DIAMOND Refseq genome short-reads alignment: $runtime min" >> \
        $OUTPUT_DIR/time/${base}_time.log
    echo ========================================== >> \
        $OUTPUT_DIR/time/${base}_time.log
fi


if $DIAMOND_SEED && [ ! -f $OUTPUT_DIR/diamond/${base}_seed.dmdout ]; then
    echo =======================================================================
    echo SEED annotation with DIAMOND ...

    if [ $rna_samples ]; then
        infile=$OUTPUT_DIR/main/${base}_mRNA.fq
    else
        infile=$OUTPUT_DIR/main/${base}_fltr.fq
    fi

    start=`date +%s`
    # Align with SEED subsystem database
    $DIAMOND blastx --threads $n_threads \
        -d $REF_DIR/subsys_db \
        -q $infile \
        -o $OUTPUT_DIR/diamond/${base}_seed.dmdout \
        -t $OUTPUT_DIR/dmnd_tmp -f 6 -k 1 --sensitive

    # $DIAMOND blastx -p $n_threads -d $REF_DIR/subsys_d \
    #     -q $infile \
    #     -a $OUTPUT_DIR/diamond/${base}.Subsys \
    #     -t ./dmnd_tmp -k 1 --sensitive
    # $DIAMOND view --daa ${base}.Subsys.daa -f 6 \
    #     -o $OUTPUT_DIR/diamond/${base}_seed.dmdout

    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "DIAMOND SEED gene annotation: $runtime min" >> \
        $OUTPUT_DIR/time/${base}_time.log
    echo ========================================== >> \
        $OUTPUT_DIR/time/${base}_time.log
fi

if $DIAMOND_COUNT; then
    echo =======================================================================
    echo Aggregate DIAMOND results and count reads...

    start=`date +%s`
    if [ ! -f $OUTPUT_DIR/counts/dmnd_RefSeq/${base}_refseq.dm_function.tsv ];
    then 
        echo Diamond RefSeq Results aggreggation ...
        python $PYSCRIPT_DIR/DIAMOND_analysis_counter.py -O \
            -I $OUTPUT_DIR/diamond/${base}_refseq.dmdout \
            -D $REF_DIR/RefSeq_bac.fa
        python $PYSCRIPT_DIR/DIAMOND_analysis_counter.py -F \
            -I $OUTPUT_DIR/diamond/${base}_refseq.dmdout \
            -D $REF_DIR/RefSeq_bac.fa
        mkdir -p $OUTPUT_DIR/counts/dmnd_RefSeq/
        mv $OUTPUT_DIR/diamond/${base}_refseq.dm_organism.tsv $OUTPUT_DIR/counts/dmnd_RefSeq/
        mv $OUTPUT_DIR/diamond/${base}_refseq.dm_function.tsv $OUTPUT_DIR/counts/dmnd_RefSeq/
        echo Completed counting of the RefSeq anotated reads!
    fi
    
    if [ ! -f $OUTPUT_DIR/counts/dmnd_NR/${base}_nr_contigs.dm_organism.tsv ];
    then 
        echo Diamond NR protein results aggreggation ...
        mkdir -p $OUTPUT_DIR/counts/dmnd_NR/
        for file in $OUTPUT_DIR/diamond/${base}_nr_*
        do 
            echo $file
            if [[ ! -s $file ]]; then
                 echo $file is empty.
                 continue
            fi
            python $PYSCRIPT_DIR/DIAMOND_analysis_counter.py -F \
                -I $file -D $REF_DIR/nr
            mv ${file%.dmdout}.dm_function.tsv $OUTPUT_DIR/counts/dmnd_NR/
        done
        echo Completed counting of the RefSeq anotated reads!
    fi

    if [ ! -f $OUTPUT_DIR/counts/dmnd_SEED/${base}_seed.receipt ];
    then
        echo Diamond SEED results aggregation ...
        mkdir -p $OUTPUT_DIR/counts/dmnd_SEED/
        python $SAMSA2/DIAMOND_subsystems_analysis_counter.py \
            -I $OUTPUT_DIR/diamond/${base}_seed.dmdout \
            -D $REF_DIR/subsys_db.fa \
            -O $OUTPUT_DIR/counts/dmnd_SEED/${base}_seed.hierarchy \
            -P $OUTPUT_DIR/counts/dmnd_SEED/${base}_seed.receipt
        # This quick program reduces down identical hierarchy annotations
        python $SAMSA2/subsys_reducer.py \
            -I $OUTPUT_DIR/counts/dmnd_SEED/${base}_seed.hierarchy
        echo Completed counting of the SEED  anotated reads!
    fi
    
    end=`date +%s`
    runtime=$(((end-start)/60))
    echo "SAMSA2 python scripts summarise DIAMOND results: $runtime min" >> \
        $OUTPUT_DIR/time/${base}_time.log
    echo ========================================== >> \
        $OUTPUT_DIR/time/${base}_time.log

fi



end0=`date +%s`
runtime0=$(((end0-start0)/60))
echo ========================================== >> $OUTPUT_DIR/time/${base}_time.log
echo "    " >> $OUTPUT_DIR/time/${base}_time.log
echo "ENTIRE WORKFLOW RUNTIME: $runtime0 min" >> \
    $OUTPUT_DIR/time/${base}_time.log

echo "  "
echo =======================================================================
echo "  "
echo Completed the workflow!
