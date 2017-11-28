#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Prepare species and gene data sets for statistical analysis, using different
## bioinformatics tools.
##
## author: sankaran.kris@gmail.com
## date: 11/27/2017

###############################################################################
## Some analysis using MIDAS
###############################################################################
library("stringr")
start_ix <- 1
end_ix <- 5

## setup paths
midas_path <- "/scratch/users/kriss1/applications/MIDAS"
python_path <- sprintf("%s:%s", Sys.getenv("PYTHONPATH"), midas_path)
path <- sprintf("%s:%s", Sys.getenv("PATH"), file.path(midas_path, "scripts"))
Sys.setenv("PYTHONPATH" = python_path)
Sys.setenv("PATH" = path)
Sys.setenv("MIDAS_DB" = file.path(midas_path, "database", "midas_db_v1.2"))

## load required modules
system("module load biology")
system("module load samtools/1.6")

outdir <- "../data/processed/"
dir.create("../data/processed/")

input_files <- list.files("../data/", "*.fq", full.names = TRUE)
input_files <- gsub("_1P.fq||_2P.fq", "", input_files)
input_files <- input_files[start_ix:end_ix]

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
