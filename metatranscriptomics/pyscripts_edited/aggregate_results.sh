BASE_DIR=$SCRATCH/Projects/perturbation_16s/
PYSCRIPT_DIR=$BASE_DIR/metatranscriptomics/pyscripts_edited
RES_DIR=${1:-$PI_SCRATCH/resilience/metatranscriptomics/processed/}
OUT=$RES_DIR/final_results/

mkdir -p $OUT
cd $RES_DIR
for subj_dir in *; do
    echo SUBJECT: $subj_dir
    echo Compute NR functions for contigs ...
    python $PYSCRIPT_DIR/generate_count_matrix.py \
        -D $RES_DIR/$subj_dir/counts/dmnd_NR/ \
        -O $OUT/${subj_dir}_contigs_nr_function.csv \
        -S _nr_contigs.dm_function.tsv
    echo Compute NR functions for unassembled ...
    python $PYSCRIPT_DIR/generate_count_matrix.py \
        -D $RES_DIR/$subj_dir/counts/dmnd_NR/ \
        -O $OUT/${subj_dir}_unassembled_nr_function.csv \
        -S _nr_unassembled.dm_function.tsv
    echo Compute RefSeq organism matrix ...
    python $PYSCRIPT_DIR/generate_count_matrix.py \
        -D $RES_DIR/$subj_dir/counts/dmnd_RefSeq/ \
        -O $OUT/${subj_dir}_refseq_organism.csv \
        -S _refseq.dm_organism.tsv
    echo Compute RefSeq function matrix ...
    python $PYSCRIPT_DIR/generate_count_matrix.py \
        -D $RES_DIR/$subj_dir/counts/dmnd_RefSeq/ \
        -O $OUT/${subj_dir}_refseq_function.csv \
        -S _refseq.dm_function.tsv
    echo Compute SEED hierarchical functions ...
    python $PYSCRIPT_DIR/generate_count_matrix.py \
        -D $RES_DIR/$subj_dir/counts/dmnd_SEED/ \
        -O $OUT/counts/${subj_dir}_SEED.csv \
        -S _seed.hierarchy.reduced
done




