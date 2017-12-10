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
parser <- add_argument(parser, "--subdir", help = "The subdirectory of data/ containing all the raw data", default = "metagenomic")
argv <- parse_args(parser)

###############################################################################
## Setup paths and load modules
###############################################################################
midas_path <- "/scratch/users/kriss1/applications/MIDAS"
python_path <- sprintf("%s:%s", Sys.getenv("PYTHONPATH"), midas_path)
path <- sprintf("%s:%s", Sys.getenv("PATH"), file.path(midas_path, "scripts"))
Sys.setenv("PYTHONPATH" = python_path)
Sys.setenv("PATH" = path)
Sys.setenv("MIDAS_DB" = file.path(midas_path, "database", "midas_db_v1.2"))

system("module load biology; module load samtools/1.6")

###############################################################################
## Define input and output directories
###############################################################################
indir <- file.path("..", "data", argv$subdir)
outdir <- file.path(subdir, "processed")
dir.create(outdir)

input_files <- list.files(indir, "*.fq", full.names = TRUE)
input_files <- unique(gsub("_1P.fq||_2P.fq", "", input_files))

## Identify and filter away procesed samples
processed_files <- list.files(outdir, full.names = TRUE)
for (f in processed_files) {
  snps_dir <- file.path(f, "snps", "output")
  if (length(list.files(snps_dir)) > 0) {
    input_files <- input_files[!grepl(basename(f), input_files)]
  }
}

input_files <- input_files[argv$start_ix:argv$end_ix]

###############################################################################
## Loop over input, performing profiling one file at a time
###############################################################################

for (f in input_files) {
  f1 <- paste0(f, "_1P.fq")
  f2 <- paste0(f, "_2P.fq")
  meas <- str_extract(f1, "M[0-9]+")

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

  ## snps profiling
  cmd <- sprintf(
    "run_midas.py snps %s/%s -1 %s -2 %s",
    outdir, meas, f1, f2
  )
  system(cmd)
}
