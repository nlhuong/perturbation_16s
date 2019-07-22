export MEM=4G
export NCPU=16

export BASE_DIR=$SCRATCH/Projects/perturbation_16s
export PI_BASE_DIR=$PI_SCRATCH/resilience/metatranscriptomics
export PYSCRIPT_DIR=$BASE_DIR/metatranscriptomics/pyscripts_edited/
export IN_DIR=${1:-$PI_BASE_DIR/processed/final_results/batches}
export OUT_DIR=${2:-$PI_BASE_DIR/processed/final_results}
export LOG_DIR=${3:-$PI_BASE_DIR/logs/final_results}
export JOB_NAME=mege-batches

mkdir -p $OUT_DIR
mkdir -p $LOG_DIR

sbatch <<EOT
#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=$JOB_NAME
#################
#a file for job output, you can check job progress, append the job ID with %j to make it unique
#SBATCH --output=${LOG_DIR}/${JOB_NAME}.%j.out
#################
# a file for errors from the job
#SBATCH --error=${LOG_DIR}/${JOB_NAME}.%j.err
#################
#time you think you need; default is 2 hours
#format could be dd-hh:mm:ss, hh:mm:ss, mm:ss, or mm
#SBATCH --time=05:00:00
#################
#Quality of Service (QOS); think of it as job priority, there is also --qos=long for with a max job length of 7 days, qos normal is 48 hours.
# REMOVE "normal" and set to "long" if you want your job to run longer than 48 hours,
# NOTE- in the hns partition the default max run time is 7 days , so you wont need to include qos

##SBATCH --qos=normal

# We are submitting to the dev partition, there are several on sherlock: normal, gpu, owners, hns, bigmem (jobs requiring >64Gigs RAM)
# I DON"T HAVE ACCESS to hns ???
#SBATCH -p normal,hns,stat,owners
#################
#number of nodes you are requesting, the more you ask for the longer you wait
##SBATCH --nodes=1
#################
#number of nodes you are requesting, the more you ask for the longer you wait
#SBATCH --mem-per-cpu=$MEM
#SBATCH --cpus-per-task=$NCPU
#################

# Have SLURM send you an email when the job ends or fails, careful, the email could end up in your clutter folder
# Also, if you submit hundreds of jobs at once you will get hundreds of emails.

#SBATCH --mail-type=END,FAIL # notifications for job done & fail
# Remember to change this to your email
#SBATCH --mail-user=$USER@stanford.edu
#################

# now run normal batch commands

echo Merging batches...

srun python $PYSCRIPT_DIR/merge_batches.py \
    $IN_DIR -outdir $OUT_DIR 
EOT
