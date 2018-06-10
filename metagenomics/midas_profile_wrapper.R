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

input_files <- list.files(argv$subdir, "*.fastq*", full.names = TRUE, recursive = TRUE)
input_files <- unique(gsub("_R1_001||_R2_001", "", input_files))
# files should start with 'M' for measure id
base_names <- sapply(input_files, function(f) basename(f))
input_files <- input_files[startsWith(base_names, "M")]
cat("Number of measures input files", length(input_files), "\n")

## Identify and filter away procesed samples
processed_files <- list.files(file.path(argv$subdir, "..", "processed"), full.names = TRUE)
cat("Number of processed files", length(processed_files), "\n")
non_empty_processed <- c()
for (f in processed_files) {
  genes_dir <- file.path(f, "genes", "output")
  if (length(list.files(genes_dir)) > 0) {
    non_empty_processed <- c(non_empty_processed, basename(f))
    input_files <- input_files[!grepl(basename(f), input_files)]
  }
}
print(length(non_empty_processed))

print("Running MIDAS for files:")
base_name <- sapply(input_files, function(f) basename(f))
names(base_name) <- NULL
print(base_name)
n_files <- length(input_files)
n_per_batch <- 4

if (n_files > 0) {
  for (i in seq(1, n_files, n_per_batch)) {
    cmd <- sprintf("bash submit.sh %s %s %s", i, min(n_files, i + n_per_batch - 1), argv$subdir)
    message(cmd)
    system(cmd)
  }
}
