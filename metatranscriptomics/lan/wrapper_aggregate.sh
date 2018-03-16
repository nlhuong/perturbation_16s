RES_DIR=${1:-$PI_SCRATCH/resilience/metatranscriptomics/processed/}
OUT=$RES_DIR/final_results/
CODE_DIR=$SCRATCH/Projects/perturbation_16s/metatranscriptomics/lan/
PARTITION=normal
mkdir -p $OUT
cd $RES_DIR
for subj_dir in *; do
    bash $CODE_DIR/aggregate_counts.sh $RES_DIR/$subj_dir $OUT $PARTITION
    sleep 1
done




