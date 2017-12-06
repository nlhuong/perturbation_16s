#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Wrapper for prepare.R, that allows samples to be processed in parallel on
## sherlock.
##
## author: sankaran.kris@gmail.com
## date: 11/27/2017

input_files <- list.files("../data/", "*.fq", full.names = TRUE)
input_files <- unique(gsub("_1P.fq||_2P.fq", "", input_files))

## Identify and filter away procesed samples
processed_files <- list.files("../data/processed", full.names = TRUE)
for (f in processed_files) {
  snps_dir <- file.path(f, "snps", "output")
  if (length(list.files(snps_dir)) > 0) {
    input_files <- input_files[!grepl(basename(f), input_files)]
  }
}

n_files <- length(input_files)
n_per_batch <- 3

for (i in seq(1, n_files, n_per_batch)) {
  cmd <- sprintf("bash submit.sh %s %s", i, min(n_files, i + n_per_batch - 1))
  message(cmd)
  system(cmd)
}
