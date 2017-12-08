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
library("feather")

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

coverage <- read_tsv("../data/merged/coverage.tsv")
depths <- read_feather("../data/merged/depths.feather")
species <- read_tsv("../data/merged/species_prevalence.tsv")
copy_num <- read_feather("../data/merged/copy_num.feather")
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
## Gene depths
###############################################################################
depths_df <- depths[, -c(1, 2)] %>%
  as.data.frame()
rownames(depths_df) <- depths$gene_id
pheatmap(t(depths_df), show_colnames = FALSE)

## the only genes that get picked up are associated with Akkermansis muciniphila
depths$species %>% unique()

###############################################################################
## Species summary statistics
###############################################################################
species %>%
  arrange(desc(mean_coverage)) %>%
  .[["species_id"]]
