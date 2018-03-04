BASE_DIR=$SCRATCH/Projects/perturbation_16s/
SUBJECT=DBUr_Sub
IN=${1:-$BASE_DIR/data/metatranscriptomics/resilience/input/$SUBJECT/}
OUT=${2:-$PI_SCRATCH/resilience/metatranscriptomics/processed/$SUBJECT/}
TIME=${3:-24:00:00}
NCPUS=${4:-20}
MEM=${5:-4G}
LOG_DIR=${5:-$BASE_DIR/logs/$SUBJECT/}

cd $BASE_DIR/metatranscriptomics/lan/
for file in $IN/*_1P.fq.gz; do
    FWD=$file
    REV=${FWD%_1P.fq.gz}_2P.fq.gz
    bash ./test.sh workflow $IN $OUT $FWD $REV $TIME $NCPUS $MEM
    sleep 1 
done
