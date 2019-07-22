CODE_DIR=$SCRATCH/Projects/perturbation_16s/metatranscriptomics/
PI_BASE_DIR=$PI_SCRATCH/resilience/metatranscriptomics
PROC_DIR=${1:-$PI_BASE_DIR/processed/}
OUT_DIR=${2:-$PI_BASE_DIR/processed/final_results/batches}

mkdir -p $OUT_DIR

cd $PROC_DIR
declare -a arr=("Relman_RNAseq_16" "Relman_RNAseq_17" "Relman_RNAseq_18" "Relman_RNAseq_21" "Arati_R_plate_11" "Arati_R_plate_9")

for subdir in "${arr[@]}"; do
    bash $CODE_DIR/aggregate_counts.sh $PROC_DIR $OUT_DIR $subdir
    sleep 1
done


