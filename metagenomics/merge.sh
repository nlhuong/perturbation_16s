export PI_BASE_DIR=$PI_SCRATCH/resilience/metagenomics
export LOG_DIR=$PI_BASE_DIR/logs
mkdir -p $LOG_DIR
#!/bin/bash
sbatch <<EOT
#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=midas_merging
#################
#a file for job output, you can check job progress, append the job ID with %j to make it unique
#SBATCH --output=$LOG_DIR/merge%j.out
#################
# a file for errors from the job
#SBATCH --error=$LOG_DIR/merge%j.err
#################
#time you think you need; default is 2 hours
#format could be dd-hh:mm:ss, hh:mm:ss, mm:ss, or mm
#SBATCH --time=15:00:00
#################
#Quality of Service (QOS); think of it as job priority, there is also --qos=long for with a max job length of 7 days, qos normal is 48 hours.
# REMOVE "normal" and set to "long" if you want your job to run longer than 48 hours,
# NOTE- in the hns partition the default max run time is 7 days , so you wont need to include qos

#SBATCH --qos=normal

# We are submitting to the dev partition, there are several on sherlock: normal, gpu, owners, hns, bigmem (jobs requiring >64Gigs RAM)
#
#SBATCH -p normal,hns
#################
#number of nodes you are requesting, the more you ask for the longer you wait
#SBATCH --mem-per-cpu=20G
#SBATCH --nodes=1
#################

# Have SLURM send you an email when the job ends or fails, careful, the email could end up in your clutter folder
# Also, if you submit hundreds of jobs at once you will get hundreds of emails.

#SBATCH --mail-type=END,FAIL # notifications for job done & fail
# Remember to change this to your email
#SBATCH --mail-user=$USER@stanford.edu

module load R/3.4.0

#now run normal batch commands
cd $SCRATCH/Projects/perturbation_16s/metagenomics
Rscript midas_merge.R
EOT
