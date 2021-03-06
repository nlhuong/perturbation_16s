#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Prepare species and gene data sets for statistical analysis, using different
## bioinformatics tools.
##
## For details about MIDAS, refer to
## https://github.com/snayfach/MIDAS/blob/master/docs/species.md
## https://github.com/snayfach/MIDAS/blob/master/docs/cnvs.md
##
## author: sankaran.kris@gmail.com
## date: 11/27/2017

library("stringr")
library("argparser")
parser <- arg_parser("Apply MIDAS profiling to raw reads")
parser <- add_argument(parser, "--start_ix", help = "Start index of files for input", default = 1)
parser <- add_argument(parser, "--end_ix", help = "End index of files for input", default = 5)
parser <- add_argument(parser, "--indir", help = "The path to the directory containing all the raw data", default = "$PI_SCRATCH/resilience/metagenomics/raw/")
argv <- parse_args(parser)

###############################################################################
## Setup paths and load modules
###############################################################################
midas_path <- file.path(Sys.getenv("SCRATCH"), "applications/MIDAS")
python_path <- sprintf("%s:%s", Sys.getenv("PYTHONPATH"), midas_path)
path <- sprintf("%s:%s", Sys.getenv("PATH"), file.path(midas_path, "scripts"))
Sys.setenv("PYTHONPATH" = python_path)
Sys.setenv("PATH" = path)
Sys.setenv("MIDAS_DB" = file.path(midas_path, "database", "midas_db_v1.2"))

###############################################################################
## Define input and output directories
###############################################################################
outdir <- file.path(argv$indir, "..", "processed")
dir.create(outdir)

input_files <- list.files(argv$indir, "*.fastq*", recursive = TRUE, full.names = TRUE)
input_files <- unique(gsub("_R1_001||_R2_001", "", input_files))
# files should start with 'M' for measure id
base_names <- sapply(input_files, function(f) basename(f))
input_files <- input_files[startsWith(base_names, "M")]

## Identify and filter away procesed samples
processed_files <- list.files(outdir, full.names = TRUE)
for (f in processed_files) {
  genes_dir <- file.path(f, "genes", "output")
  if (length(list.files(genes_dir)) > 0) {
    input_files <- input_files[!grepl(basename(f), input_files)]
  }
}

input_files <- input_files[argv$start_ix:argv$end_ix]
print(input_files)

###############################################################################
## Loop over input, performing profiling one file at a time
###############################################################################
for (f in input_files) {
  f1 <- gsub(".fastq", "_R1_001.fastq", f)
  f2 <- gsub(".fastq", "_R2_001.fastq", f)
  meas <- str_extract(f1, "M[0-9]+")
  if (any(is.na(meas), meas == "NA")){
    #meas <- f
    next 
  }

  ## species profiling
  cmd <- sprintf(
    "run_midas.py species %s/%s -1 %s -2 %s",
    outdir, meas, f1, f2
  )
  system(cmd)

  ## gene profiling
  cmd <- sprintf(
    "run_midas.py genes %s/%s -1 %s -2 %s",
    outdir, meas, f1, f2
  )
  system(cmd)
}
