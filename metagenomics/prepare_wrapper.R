#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Wrapper for prepare.R, that allows samples to be processed in parallel on
## sherlock.
##
## author: sankaran.kris@gmail.com
## date: 11/27/2017

input_files <- list.files("../data/", "*.fq", full.names = TRUE)
input_files <- gsub("_1P.fq||_2P.fq", "", input_files)
## n_files <- length(input_files)
n_files <- 15
n_per_batch <- 5

start_ix <- seq(1, n_files, n_per_batch)

for (i in start_ix) {
  for (j in i:(i + n_per_batch - 1)) {
    system(sprintf("bash submit.sh %s %s", i, j))
  }
}
