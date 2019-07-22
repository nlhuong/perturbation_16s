export JOB_NAME=resilience_fpca
export NCPUS=16
export MEM=8GB
sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=$JOB_NAME
##SBATCH --begin=now+7days
#SBATCH --dependency=singleton
#SBATCH --time=03:00:00
#SBATCH -p normal,hns,stat,owners

#SBATCH --mail-type=FAIL
#SBATCH --output=${JOB_NAME}.%j.out
#SBATCH --error=${JOB_NAME}.%j.err
#SBATCH --mem-per-cpu=$MEM
#SBATCH --cpus-per-task=$NCPUS

## Insert the command to run below.

Rscript fclust.R

## Resubmit the job for the next execution
#sbatch $0

EOT

