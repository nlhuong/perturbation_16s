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
export base_dir=/scratch/users/lanhuong/Projects/perturbation_16s/
export data_dir=$base_dir/data/metatranscriptomics/resilience
export processed_dir=$data_dir/output
export raw_dir=$data_dir/input

cd $base_dir/metatranscriptomics/pyscripts_edited

srun -p normal -c 10 python summarize_stats.py \
    -R $raw_dir -P $processed_dir -O "./DBUr_stats.csv" -N 10

