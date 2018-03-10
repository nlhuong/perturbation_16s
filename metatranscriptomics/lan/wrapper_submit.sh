SUBJECT=EAYr_Sub
BASE_DIR=$SCRATCH/Projects/perturbation_16s/
SCRIPT=$BASE_DIR/metatranscriptomics/lan/submit.sh
IN=${1:-$BASE_DIR/data/metatranscriptomics/resilience/input/$SUBJECT/}
OUT=${2:-$PI_SCRATCH/resilience/metatranscriptomics/processed/$SUBJECT/}
TIME=${3:-20:00:00}
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
