export NCPU=24
export filename=all_samples

export BASE_DIR=$SCRATCH/Projects/perturbation_16s
export PI_BASE_DIR=$PI_SCRATCH/resilience/metagenomics
export RAW_DIR=$PI_BASE_DIR/raw
export MERGED_FILE=$PI_BASE_DIR/merged/count_reads.txt
export CODE_DIR=$BASE_DIR/metagenomics
export LOG_DIR=$PI_BASE_DIR/logs/stats
export OUT_DIR=$PI_BASE_DIR/stats
export JOB_NAME=midas-stats

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
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=$USER@stanford.edu

module load python

srun python $CODE_DIR/midas_summary_stats.py \
    $RAW_DIR $MERGED_FILE -outfile ${OUT_DIR}/${filename}_stats.csv \
    -ncores $NCPU 

EOT
