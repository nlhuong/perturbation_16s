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

#SUBDIR=Arati_R_plate_11
SUBDIR=Relman_RNAseq_21

IN=${1:-$PI_BASE_DIR/raw/$SUBDIR}
OUT=${2:-$PI_BASE_DIR/processed/$SUBDIR}
LOG=${3:-$PI_BASE_DIR/logs/$SUBDIR/sh_uniref/sh_count}

TIME=${3:-24:00:00}
NCPUS=${4:-8}
MEM=${5:-8G}


cd $IN
mkdir -p $OUT
running_jobs=($(squeue -u $USER | awk '{ gsub("wf-","",$3); print $3 }'))
for fwd_file in *_R1_001.fastq; do
    base="$(echo $fwd_file | cut -d '_' -f1-3)"
    base_running=false
    for job in "${running_jobs[@]}"; do
        if [[ "$base" == *${job}* ]]; then
            base_running=true
        fi
    done
    if [[ $base_running == true ]]; then
        echo "sample ${base} already running." 
        continue
    fi
    if [ ! -s $OUT/counts/dmnd_UniRef/${base}_uniref_abund.csv ]; then
        echo $base
        FWD=$fwd_file
        REV=${FWD%_R1_001.fastq}_R2_001.fastq
        bash $CODE_DIR/submit.sh $IN $OUT $FWD $REV $LOG $TIME $NCPUS $MEM
        sleep 1 
    fi
done

#done
