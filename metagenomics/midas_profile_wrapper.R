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
parser <- add_argument(parser, "--subdir", help = "The subdirectory of data/ containing all the raw data", default = "metagenomic")

indir <- file.path("..", "data", argv$subdir)
input_files <- list.files(indir, "*.fq", full.names = TRUE)
input_files <- unique(gsub("_1P.fq||_2P.fq", "", input_files))

## Identify and filter away procesed samples
processed_files <- list.files(file.path(indir, "processed"), full.names = TRUE)
for (f in processed_files) {
  snps_dir <- file.path(f, "snps", "output")
  if (length(list.files(snps_dir)) > 0) {
    input_files <- input_files[!grepl(basename(f), input_files)]
  }
}

n_files <- length(input_files)
n_per_batch <- 4

for (i in seq(1, n_files, n_per_batch)) {
  cmd <- sprintf("bash submit.sh %s %s %s", i, min(n_files, i + n_per_batch - 1), indir)
  message(cmd)
  system(cmd)
}
