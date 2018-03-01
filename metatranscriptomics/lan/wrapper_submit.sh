BASE_DIR=$SCRATCH/Projects/perturbation_16s/

IN=${1:-$BASE_DIR/data/metatranscriptomics/resilience/input/DBUr_Sub/}
OUT=${2:-$BASE_DIR/data/metatranscriptomics/resilience/output/DBUr_Sub/}
TIME=${3:-24:00:00}
NCPUS=${4:-10}
MEM=${5:-4G}

cd $BASE_DIR/metatranscriptomics/lan/
for file in $IN/*_1P.fq.gz; do
    FWD=$file
    REV=${FWD%_1P.fq.gz}_2P.fq.gz
    bash ./test.sh workflow $IN $OUT $FWD $REV $TIME $NCPUS $MEM
    sleep 1 
done
