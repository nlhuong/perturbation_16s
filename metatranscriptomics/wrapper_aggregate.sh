CODE_DIR=$SCRATCH/Projects/perturbation_16s/metatranscriptomics/
PI_BASE_DIR=$PI_SCRATCH/resilience/metatranscriptomics
RES_DIR=${1:-$PI_BASE_DIR/processed/}
OUT=$PI_BASE_DIR/logs/final_results/
mkdir -p $OUT
cd $RES_DIR
for subj_dir in *; do
    bash $CODE_DIR/aggregate_counts.sh $RES_DIR/$subj_dir $OUT
    sleep 1
done




