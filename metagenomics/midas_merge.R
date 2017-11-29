#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Merge profiling results from midas into single count tables for further
## statistical analysis.
##
## author: sankaran.kris@gmail.com
## date: 11/29/2017

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
samples_dir <- file.path("..", "data", "processed")
merged_dir <- file.path("..", "data", "merged")
genes_dirs <- file.path(merged_dir, "genes")
dir.create(merged_dir)
dir.create(genes_dirs)

base_cmd <- "merge_midas.py %s %s -i %s -t dir"
system(sprintf(base_cmd, "genes", merged_dir, samples_dir))
system(sprintf(base_cmd, "genes", genes_dirs, samples_dir))
system(sprintf(base_cmd, "snps", genes_dirs, samples_dir))

###############################################################################
## Combine gene coverage and depth data
###############################################################################
