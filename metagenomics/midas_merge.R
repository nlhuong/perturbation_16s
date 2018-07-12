#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Merge profiling results from midas into single count tables for further
## statistical analysis.
##
## See https://github.com/snayfach/MIDAS/blob/master/docs/merge_cnvs.md for
## underlying script details.
##
## author: sankaran.kris@gmail.com
## date: 11/29/2017

library("tidyverse")
library("feather")
library("argparser")

parser <- arg_parser("Merge MIDAS species and gene results across samples")
parser <- add_argument(parser, "--subdir", help = "The directory containing the processed data", default = file.path(Sys.getenv("PI_SCRATCH"), "resilience", "metagenomics"))
argv <- parse_args(parser)

bind_wrapper <- function(x) {
  x %>%
    bind_rows %>%
    select(species, gene_id, starts_with("M"))
}

###############################################################################
## Setup PATH to MIDAS
###############################################################################
midas_path <- file.path(Sys.getenv("SCRATCH"), "applications/MIDAS")
python_path <- sprintf("%s:%s", Sys.getenv("PYTHONPATH"), midas_path)
path <- sprintf("%s:%s", Sys.getenv("PATH"), file.path(midas_path, "scripts"))
Sys.setenv("PYTHONPATH" = python_path)
Sys.setenv("PATH" = path)
Sys.setenv("MIDAS_DB" = file.path(midas_path, "database", "midas_db_v1.2"))

###############################################################################
## Merge data across samples
###############################################################################
samples_dir <- file.path(argv$subdir, "processed")
merged_dir <- file.path(argv$subdir, "merged")
species_dir <- file.path(merged_dir, "species")
genes_dirs <- file.path(merged_dir, "genes")
dir.create(merged_dir)
dir.create(species_dir)
dir.create(genes_dirs)

base_cmd <- "merge_midas.py %s %s -i %s -t dir"
system(sprintf(base_cmd, "species", species_dir, samples_dir))
system(sprintf(base_cmd, "genes", genes_dirs, samples_dir))

###############################################################################
## Combine gene coverage and depth data
###############################################################################
genes_f <- list.files(genes_dirs, full.names = TRUE)
reads <- list()
depths <- list()
copy_num <- list()

for (i in seq_along(genes_f)) {
  message("Merging ", genes_f[i])

  reads[[i]] <- read_tsv(file.path(genes_f[i], "genes_reads.txt")) 
  reads[[i]]$species <- basename(genes_f[i])

  depths[[i]] <- read_tsv(file.path(genes_f[i], "genes_depth.txt"))
  depths[[i]]$species <- basename(genes_f[i])

  copy_num[[i]] <- read_tsv(file.path(genes_f[i], "genes_copynum.txt"))
  copy_num[[i]]$species <- basename(genes_f[i])
}

###############################################################################
## Write results to file
###############################################################################
write_feather(bind_wrapper(reads), file.path(merged_dir, "gene_reads.feather"))
write_feather(bind_wrapper(depths), file.path(merged_dir, "gene_depths.feather"))
write_feather(bind_wrapper(copy_num), file.path(merged_dir, "gene_copy_num.feather"))
file.copy(file.path(species_dir, "coverage.txt"), 
    file.path(merged_dir, "species_identifier_coverage.tsv"))
