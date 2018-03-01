export BASE_DIR=$SCRATCH/Projects/perturbation_16s
export SCRIPT=${1:-workflow}
export IN=${2:-$BASE_DIR/data/metatranscriptomics/resilience/input/DBUr_Sub/}
export OUT=${3:-$BASE_DIR/data/metatranscriptomics/resilience/output/DBUr_Sub/}
export FWD=${4:-M3303_DBUsw_2r_TrM31_1P.fq.gz}
export REV=${5:-M3303_DBUsw_2r_TrM31_2P.fq.gz}
export TIME=${6:-24:00:00}
export NCPUS=${7:-10}
export MEM=${8:-4G}
sbatch <<EOT
#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=$SCRIPT
#################
#a file for job output, you can check job progress, append the job ID with %j to make it unique
#SBATCH --output=../../logs/${SCRIPT}.out
#################
# a file for errors from the job
#SBATCH --error=../../logs/${SCRIPT}.err
#################
#time you think you need; default is 2 hours
#format could be dd-hh:mm:ss, hh:mm:ss, mm:ss, or mm
#SBATCH --time=$TIME
#################
#Quality of Service (QOS); think of it as job priority, there is also --qos=long for with a max job length of 7 days, qos normal is 48 hours.
# REMOVE "normal" and set to "long" if you want your job to run longer than 48 hours,
# NOTE- in the hns partition the default max run time is 7 days , so you wont need to include qos

#SBATCH --qos=normal

# We are submitting to the dev partition, there are several on sherlock: normal, gpu, owners, hns, bigmem (jobs requiring >64Gigs RAM)
# I DON"T HAVE ACCESS to hns ???
#SBATCH -p normal
#################
#number of nodes you are requesting, the more you ask for the longer you wait
#SBATCH --nodes=1
#################
#number of nodes you are requesting, the more you ask for the longer you wait
#SBATCH --mem-per-cpu=$MEM
#SBATCH --cpus-per-task=$NCPUS
#################

# Have SLURM send you an email when the job ends or fails, careful, the email could end up in your clutter folder
# Also, if you submit hundreds of jobs at once you will get hundreds of emails.

#SBATCH --mail-type=END,FAIL # notifications for job done & fail
# Remember to change this to your email
#SBATCH --mail-user=$USER@stanford.edu
#################

# now run normal batch commands
echo JOB PARAMETERS
echo =======================================================================
echo TIME REQUESTED: $TIME
echo NUMBER OF CPUS: $NCPUS
echo MEMORY PER CPU: $MEM
echo INPUT DIR: $IN
echo OUTPUT DIR: $OUT
echo " "
echo =======================================================================
echo " "

echo SAMPLE: ${FWD%_1P.fq.gz}
echo FWD FILE         = $FWD
echo REV FILE         = $REV

EOT
