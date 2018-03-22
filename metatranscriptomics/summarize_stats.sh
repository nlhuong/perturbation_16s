export BASE_DIR=$PI_SCRATCH/resilience/metatranscriptomics/
export CODE_DIR=$SCRATCH/Projects/perturbation_16s/metatranscriptomics/
export PYSCRIPT=$CODE_DIR/pyscripts_edited/summarize_stats.py

export RAW_DIR=${1:-$BASE_DIR/raw}
export PROCESSED_DIR=${2:-$BASE_DIR/processed}
export OUTPUT=${3:-$PROCESSED_DIR/metatrans_stats.csv}
export LOG_DIR=${4:-$BASE_DIR/logs}
mkdir -p $LOG_DIR

export TIME=${5:-05:00:00}
export NCPUS=${6:-8}
export MEM=${7:-4G}

export JOB_NAME=stats-metatrans
sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=$JOB_NAME
#SBATCH --output=${LOG_DIR}/${JOB_NAME}.%j.out
#SBATCH --error=${LOG_DIR}/${JOB_NAME}.%j.err
#format could be dd-hh:mm:ss, hh:mm:ss, mm:ss, or mm
#SBATCH --time=$TIME
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=$MEM
#SBATCH --cpus-per-task=$NCPUS

module load python

cd $CODE_DIR/pyscripts_edited
srun python summarize_stats.py \
    -R $RAW_DIR -P $PROCESSED_DIR -O $OUTPUT -N $NCPUS

EOT
