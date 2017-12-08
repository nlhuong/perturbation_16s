#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Merge profiling results from midas into single count tables for further
## statistical analysis.
##
## author: sankaran.kris@gmail.com
## date: 11/29/2017

library("tidyverse")
library("feather")

bind_wrapper <- function(x) {
  x %>%
    bind_rows %>%
    select(species, gene_id, starts_with("M"))
}

###############################################################################
## Setup PATH to MIDAS
###############################################################################
midas_path <- "/scratch/users/kriss1/applications/MIDAS"
python_path <- sprintf("%s:%s", Sys.getenv("PYTHONPATH"), midas_path)
path <- sprintf("%s:%s", Sys.getenv("PATH"), file.path(midas_path, "scripts"))
Sys.setenv("PYTHONPATH" = python_path)
Sys.setenv("PATH" = path)
Sys.setenv("MIDAS_DB" = file.path(midas_path, "database", "midas_db_v1.2"))

###############################################################################
## Merge data across samples
###############################################################################
samples_dir <- file.path("..", "data", "metagenomic", "processed")
merged_dir <- file.path("..", "data", "metagenomic", "merged")
genes_dirs <- file.path(merged_dir, "genes")
dir.create(merged_dir)
dir.create(genes_dirs)

base_cmd <- "merge_midas.py %s %s -i %s -t dir"
system(sprintf(base_cmd, "genes", genes_dirs, samples_dir))
system(sprintf(base_cmd, "snps", genes_dirs, samples_dir))

###############################################################################
## Combine gene coverage and depth data
###############################################################################
genes_f <- list.files(genes_dirs, full.names = TRUE)
depths <- list()
copy_num <- list()
snps <- list()
snps_info <- list()

for (i in seq_along(genes_f)) {
  message("Merging ", genes_f[i])
  depths[[i]] <- read_tsv(file.path(genes_f[i], "genes_depth.txt"))
  depths[[i]]$species <- basename(genes_f[i])

  copy_num[[i]] <- read_tsv(file.path(genes_f[i], "genes_copynum.txt"))
  copy_num[[i]]$species <- basename(genes_f[i])

  snps[[i]] <- read_tsv(file.path(genes_f[i], "snps_depth.txt"))
  snps[[i]]$species <- basename(genes_f[i])

  snps_info[[i]] <- read_tsv(file.path(genes_f[i], "snps_info.txt"))
  snps_info[[i]]$species <- basename(genes_f[i])
}

###############################################################################
## Write results to file
###############################################################################
write_feather(bind_wrapper(depths), file.path(merged_dir, "depths.feather"))
write_feather(bind_wrapper(copy_num), file.path(merged_dir, "copy_num.feather"))
write_feather(bind_wrapper(snps), file.path(merged_dir, "snps.feather"))
write_feather(bind_wrapper(snps_info), file.path(merged_dir, "snps_info.feather"))
