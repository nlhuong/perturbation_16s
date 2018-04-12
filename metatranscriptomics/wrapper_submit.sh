# Arati_R_plate_11/
# Arati_R_plate_9/
# Relman_RNAseq_21/

# Relman_RNAseq_16/
# Relman_RNAseq_17/
# Relman_RNAseq_18/

CODE_DIR=$SCRATCH/Projects/perturbation_16s/metatranscriptomics
PI_BASE_DIR=$PI_SCRATCH/resilience/metatranscriptomics

#declare -a arr=("Relman_RNAseq_16" "Relman_RNAseq_17/" "Relman_RNAseq_18")
#for SUBDIR in "${arr[@]}"
#do

SUBDIR=Arati_R_plate_11
#Relman_RNAseq_21

IN=${1:-$PI_BASE_DIR/raw/$SUBDIR}
OUT=${2:-$PI_BASE_DIR/processed/$SUBDIR}
LOG=${3:-$PI_BASE_DIR/logs/$SUBDIR/sh2_first}

TIME=${3:-18:30:00}
NCPUS=${4:-8}
MEM=${5:-6G}

cd $IN
mkdir -p $OUT
for fwd_file in *_R1_001.fastq; do
    base="$(echo $fwd_file | cut -d '_' -f1-3)"
    if [ ! -s $OUT/counts/dmnd_SEED/${base}_seed.hierarchy.reduced ]; then
        echo $base
        FWD=$fwd_file
        REV=${FWD%_R1_001.fastq}_R2_001.fastq
        bash $CODE_DIR/submit.sh $IN $OUT $FWD $REV $LOG $TIME $NCPUS $MEM
        sleep 1 
    fi
done

#done
