export BASE_DIR=$SCRATCH/Projects/perturbation_16s
export PYSCRIPT_DIR=$BASE_DIR/metatranscriptomics/pyscripts_edited/
export RES_DIR=${1:-$PI_SCRATCH/resilience/metatranscriptomics/processed/DBUr_Sub/}
export OUT=${2:-$RES_DIR/../final_results/}
export PARTITION=${3:-dev}
export LOG_DIR=$BASE_DIR/logs/
export DIR=$(basename $RES_DIR)
export JOB_NAME=agg-cnts-${DIR}
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
#SBATCH --time=01:00:00
#################
#Quality of Service (QOS); think of it as job priority, there is also --qos=long for with a max job length of 7 days, qos normal is 48 hours.
# REMOVE "normal" and set to "long" if you want your job to run longer than 48 hours,
# NOTE- in the hns partition the default max run time is 7 days , so you wont need to include qos

##SBATCH --qos=normal

# We are submitting to the dev partition, there are several on sherlock: normal, gpu, owners, hns, bigmem (jobs requiring >64Gigs RAM)
# I DON"T HAVE ACCESS to hns ???
#SBATCH -p $PARTITION
#################
#number of nodes you are requesting, the more you ask for the longer you wait
##SBATCH --nodes=1
#################
#number of nodes you are requesting, the more you ask for the longer you wait
##SBATCH --mem-per-cpu=4G
##SBATCH --cpus-per-task=1
#################

# Have SLURM send you an email when the job ends or fails, careful, the email could end up in your clutter folder
# Also, if you submit hundreds of jobs at once you will get hundreds of emails.

#SBATCH --mail-type=END,FAIL # notifications for job done & fail
# Remember to change this to your email
#SBATCH --mail-user=$USER@stanford.edu
#################

# now run normal batch commands

echo DIRECTORY: $DIR
echo Compute NR functions for contigs ...
python $PYSCRIPT_DIR/generate_count_matrix.py \
    -D $RES_DIR/counts/dmnd_NR/ \
    -O $OUT/${DIR}_contigs_nr_function.csv \
    -S _nr_contigs.dm_function.tsv

echo Compute NR functions for unassembled ...
python $PYSCRIPT_DIR/generate_count_matrix.py \
    -D $RES_DIR/counts/dmnd_NR/ \
    -O $OUT/${DIR}_unassembled_nr_function.csv \
    -S _nr_unassembled.dm_function.tsv
echo Compute RefSeq organism matrix ...
python $PYSCRIPT_DIR/generate_count_matrix.py \
    -D $RES_DIR/counts/dmnd_RefSeq/ \
    -O $OUT/${DIR}_refseq_organism.csv \
    -S _refseq.dm_organism.tsv

echo Compute RefSeq function matrix ...
python $PYSCRIPT_DIR/generate_count_matrix.py \
    -D $RES_DIR/counts/dmnd_RefSeq/ \
    -O $OUT/${DIR}_refseq_function.csv \
    -S _refseq.dm_function.tsv

echo Compute SEED hierarchical functions ...
python $PYSCRIPT_DIR/generate_count_matrix.py \
    -D $RES_DIR/counts/dmnd_SEED/ \
    -O $OUT/${DIR}_SEED.csv \
    -S _seed.hierarchy.reduced

EOT
