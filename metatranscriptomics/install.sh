#!/usr/bin/env bash
##
## Installing things for the workflow in
## https://github.com/ParkinsonLab/2017-Microbiome-Workshop
##
## author: nlhuong90@gmail.com
## date: 2/18/2018

n_threads=10

## DIRECTORIES ----------
if [ -z ${SCRATCH+x} ]; then
    #the variable $SCRATCH is unset
    echo Working on ICME cluster
    module load gcc/gcc6
    module load cmake/cmake-3.7.2
    BASE_DIR=~/Projects/perturbation_16s
    APP_DIR=~/.local/bin/
    PYSCRIPT_DIR=$BASE_DIR/metatranscriptomics/pyscripts_edited
else
    echo Working on SHERLOCK cluster
    pip2.7 install --user scipy numpy matplotlib xlwt joblib pandas
    BASE_DIR=$SCRATCH/Projects/perturbation_16s
    APP_DIR=$SCRATCH/applications/bin/
    PYSCRIPT_DIR=$BASE_DIR/metatranscriptomics/pyscripts
    if [ $SHERLOCK == "2" ]; then
        ## The following modules must be preloaded:
        # module load python/2.7.13
        # module load py-biopython/1.70
        module load biology
        module load bwa/0.7.17
        module load samtools/1.6
        module load ncbi-blast+/2.6.0
        module load cmake/3.8.1
    fi
fi

cd $APP_DIR

## fastqc
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
unzip fastqc_v0.11.7.zip
rm fastqc_v0.11.7.zip
cd FastQC
chmod +x fastqc
cd $APP_DIR

## trimmomatic
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
unzip Trimmomatic-0.36.zip
rm Trimmomatic-0.36.zip
chmod -R +x Trimmomatic-0.36

## vsearch
wget https://github.com/torognes/vsearch/releases/download/v2.7.0/vsearch-2.7.0-linux-x86_64.tar.gz
tar xzf vsearch-2.7.0-linux-x86_64.tar.gz
rm vsearch-2.7.0-linux-x86_64.tar.gz
chmod -R +x vsearch-2.7.0-linux-x86_64

## cd-hit
git clone https://github.com/weizhongli/cdhit.git
cd cdhit
make
cd cd-hit-auxtools
make
cd $APP_DIR

## blat
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat
chmod +x blat

## Infernal
wget http://eddylab.org/infernal/infernal-1.1.2-linux-intel-gcc.tar.gz
tar -zxvf infernal-1.1.2-linux-intel-gcc.tar.gz
rm infernal-1.1.2-linux-intel-gcc.tar.gz

## kaiju
git clone https://github.com/bioinformatics-centre/kaiju.git
cd kaiju/src
make

mkdir -p $DB_DIR/kaijudb
cd $DB_DIR/kaijudb
## The following steps seem to take INSANE amount of TIME there seem to be
## and issue how makeDB.sh, mkbwt and mkfmi handle parallelization.
$APP_DIR/kaiju/bin/makeDB.sh -r -t $n_threads # make NCBI reference DB
## $APP_DIR/kaiju/bin/makeDB.sh -r -t $n_threads --noDL # if already downloaded
## Make custom library?
#$APP_DIR/kaiju/bin/mkbwt -o kaiju_db -nThreads $n_threads kaiju_db.faa
#$APP_DIR/kaiju/bin/mkfmi kaiju_db

## spades
cd $APP_DIR
wget http://cab.spbu.ru/files/release3.11.1/SPAdes-3.11.1-Linux.tar.gz
tar -zxvf SPAdes-3.11.1-Linux.tar.gz
rm SPAdes-3.11.1-Linux.tar.gz

## diamond
wget http://github.com/bbuchfink/diamond/releases/download/v0.9.18/diamond-linux64.tar.gz
tar -xzf diamond-linux64.tar.gz
rm diamond-linux64.tar.gz

## sortmerna
git clone https://github.com/biocore/sortmerna.git
cd sortmerna
mkdir -p build/Release
pushd build/Release
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release ../..
make
cd $APP_DIR

## samsa2
git clone https://github.com/transcript/samsa2.git
# NOTE THAT the python_scripts/DIAMOND_analysis_counter.py have indexing errors
# Need to change lines 148, 151, 156, 163  to have [0] and [1] indices
# as opposed to [1] and [2]


## biology --------------------------------------------------------------------------------
if [ -z ${SCRATCH+x} ] || [ $SHERLOCK == "1" ]; then
    ## bwa
    conda install -c bioconda bwa

    ## samtools
    conda install -c bioconda samtools

    ## blast
    conda install -c bioconda blast

    ## biopython
    conda install -c conda-forge biopython #1.70

    echo "export PATH=$APP_DIR/vsearch-2.7.0-linux-x86_64/bin:\$PATH" >> ~/.bashrc
    echo "export PATH=$APP_DIR/Trimmomatic-0.36:\$PATH" >> ~/.bashrc
    echo "export PATH=$APP_DIR/FastQC:\$PATH" >> ~/.bashrc
    echo "export PATH=$APP_DIR/sortmerna/build/Release/src/indexdb:$APP_DIR/sortmerna/build/Release/src/sortmerna:\$PATH" >> ~/.bashrc
    echo "export PATH=$APP_DIR:\$PATH" >> ~/.bashrc
    source ~/.bashrc
fi
