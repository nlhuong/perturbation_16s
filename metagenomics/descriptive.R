#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## A first pass at visualizing the metagenomics data produced by using MIDAS
## with the default settings.
##
## author: sankaran.kris@gmail.com
## date: 12/08/2017

###############################################################################
## Setup and data
###############################################################################
library("tidyverse")
library("readxl")

opts <- list(
  "filter" = list(
    "cov" = list(
      "k" = 0.2,
      "a" = 0
    ),
    "species" = list(
      "k" = 0.2,
      "a" = 0
    )
  )
)

coverage <- read_tsv("../data/merged/coverage.txt")
depths <- read_tsv("../data/merged/depths.txt")
species <- read_tsv("../data/merged/species_prevalence.txt")
copy_num <- read_tsv("../data/merged/copy_num.txt")
mapping <- read_xlsx("../data/Mapping_Files_7bDec2017.xlsx", skip = 2)

###############################################################################
## Study the species COG coverage data
###############################################################################
keep_ix <- rowMeans(coverage[, -1] > opts$filter$cov$k) > opts$filter$cov$a
coverage <- coverage[keep_ix, ]
coverage[, -1] <- asinh(coverage[, -1])

coverage <- as.data.frame(coverage)
rownames(coverage) <- coverage$species_id
pheatmap(t(coverage[, -1]), fontsize = 6)

###############################################################################
## Study the gene depths
###############################################################################
