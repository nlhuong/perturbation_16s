#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Wrapper for prepare.R, that allows samples to be processed in parallel on
## sherlock.
##
## author: sankaran.kris@gmail.com
## date: 11/27/2017

library("argparser")
parser <- arg_parser("Wrap metagenomic and metatranscriptomic profiling")
parser <- add_argument(parser, "--workdir", help = "The directory within which to run the R process", default = "metagenomics")
parser <- add_argument(parser, "--subdir", help = "The subdirectory containing all the raw data", default = file.path(Sys.getenv("PI_SCRATCH"), "resilience", "metagenomics", "raw"))
parser <- add_argument(parser, "--logdir", help = "Relative path to directory within which to store all cluster logs", default = file.path(Sys.getenv("PI_SCRATCH"), "resilience", "metagenomics", "logs"))
argv <- parse_args(parser)

message(argv$workdir)
setwd(argv$workdir)
dir.create(argv$logdir)

input_files <- list.files(argv$subdir, "*.fq*", full.names = TRUE, recursive = TRUE)
print(input_files)
input_files <- unique(gsub("_1P||_2P", "", input_files))

## Identify and filter away procesed samples
processed_files <- list.files(file.path(argv$subdir, "processed"), full.names = TRUE)
for (f in processed_files) {
  genes_dir <- file.path(f, "genes", "output")
  if (length(list.files(genes_dir)) > 0) {
    input_files <- input_files[!grepl(basename(f), input_files)]
  }
}

n_files <- length(input_files)
n_per_batch <- 4

for (i in seq(1, n_files, n_per_batch)) {
  cmd <- sprintf("bash submit.sh %s %s %s", i, min(n_files, i + n_per_batch - 1), argv$subdir)
  message(cmd)
  system(cmd)
}
