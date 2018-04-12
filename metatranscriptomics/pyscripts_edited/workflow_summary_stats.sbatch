#!/bin/bash
#SBATCH --job-name=summarize-stats
#SBATCH --output=summarize-stats.%j.out
#SBATCH --error=summarize-stats.%j.err
#format could be dd-hh:mm:ss, hh:mm:ss, mm:ss, or mm
#SBATCH --time=06:00:00
#SBATCH -p normal
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=10

module load python

ncpu=10
subdir=Relman_RNAseq_16

BASE_DIR=$SCRATCH/Projects/perturbation_16s
PI_BASE_DIR=$PI_SCRATCH/resilience/metatranscriptomics
CODE_DIR=$BASE_DIR/metatranscriptomics/pyscripts_edited

processed_dir=$PI_BASE_DIR/processed
raw_dir=$PI_BASE_DIR/raw

srun -p normal,owners,hns -c $NCPU python $CODE_DIR/workflow_summary_stats.py \
    $raw_dir $processed_dir -outfile ${subdir}_stats.csv \
    -ncores $ncpu -subdir $subdir
