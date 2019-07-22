export NCPU=16
#export SUBDIR=Relman_RNAseq_16
export filename=all_samples

export BASE_DIR=$SCRATCH/Projects/perturbation_16s
export PI_BASE_DIR=$PI_SCRATCH/resilience/metatranscriptomics
export CODE_DIR=$BASE_DIR/metatranscriptomics/pyscripts_edited
export PROC_DIR=$PI_BASE_DIR/processed
export RAW_DIR=$PI_BASE_DIR/raw
export LOG_DIR=$PI_BASE_DIR/logs/stats
export OUT_DIR=$PI_BASE_DIR/stats
export JOB_NAME=workflow-stats

mkdir -p $OUT_DIR
mkdir -p $LOG_DIR

sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=$JOB_NAME
#SBATCH --output=${LOG_DIR}/${JOB_NAME}.%j.out
#SBATCH --error=${LOG_DIR}/${JOB_NAME}.%j.err
#format could be dd-hh:mm:ss, hh:mm:ss, mm:ss, or mm
#SBATCH --time=06:00:00
#SBATCH -p normal,owners,hns,stat
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=$NCPU

module load python

srun python $CODE_DIR/workflow_summary_stats.py \
    $RAW_DIR $PROC_DIR -outfile ${OUT_DIR}/${filename}_stats.csv \
    -ncores $NCPU 

#-subdir $SUBDIR
EOT
