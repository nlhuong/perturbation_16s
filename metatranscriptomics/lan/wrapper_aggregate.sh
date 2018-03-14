RES_DIR=${1:-$PI_SCRATCH/resilience/metatranscriptomics/processed/}
OUT=$RES_DIR/final_results/

mkdir -p $OUT
cd $RES_DIR
for subj_dir in *; do
    bash aggregate_counts.sh $RES_DIR/$subj_dir $OUT
    sleep 1
done




