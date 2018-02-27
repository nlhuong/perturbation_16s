export icme_clust=${1:0}

if [ $icme_clust ]; then
    module load gcc/gcc6 # on ICME clusters
    export BASE_DIR=~/Projects/perturbation_16s
    export APP_DIR=~/.local/bin/
    export PYSCRIPT_DIR=$BASE_DIR/metatranscriptomics/pyscripts_edited
else 
    ## The following modules must be preloaded:
    module load python/2.7.13
    module load py-biopython/1.70
    module load biology
    module load bwa/0.7.17
    module load samtools/1.6
    module load ncbi-blast+/2.6.0
    export BASE_DIR=$SCRATCH/Projects/perturbation_16s
    export APP_DIR=$SCRATCH/applications/bin/
    export PYSCRIPT_DIR=$BASE_DIR/metatranscriptomics/pyscripts
fi

# Reference directories
export REF_DIR=$BASE_DIR/data/databases
export KAIJUBD_DIR=$REF_DIR/kaijudb

# SortMeRNA directories
export SORTMERNA_DIR=$APP_DIR/sortmerna

# Main workflow code
export CODE=$BASE_DIR/metatranscriptomics/lan
