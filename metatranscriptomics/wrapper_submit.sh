# Arati_R_plate_11/
# Arati_R_plate_9/
# Relman_RNAseq_16/
# Relman_RNAseq_17/
# Relman_RNAseq_18/
# Relman_RNAseq_21/

SUBJECT=DBUr_Sub
BASE_DIR=$SCRATCH/Projects/perturbation_16s/
SCRIPT=$BASE_DIR/metatranscriptomics/submit.sh
IN=${1:-$BASE_DIR/data/metatranscriptomics/resilience/input/$SUBJECT/}
OUT=${2:-$PI_SCRATCH/resilience/metatranscriptomics/processed/$SUBJECT/}
TIME=${3:-24:00:00}
NCPUS=${4:-8}
MEM=${5:-4G}
mkdir -p $OUT
cd $IN
for file in *_1P.fq.gz; do
    base=${file%_1P.fq.gz}
    if [ ! -f $OUT/counts/dmnd_SEED/${base}* ]; then
        echo $base
        FWD=$file
        REV=${FWD%_1P.fq.gz}_2P.fq.gz
        bash $SCRIPT workflow $IN $OUT $FWD $REV $TIME $NCPUS $MEM
        sleep 1 
    fi
done
