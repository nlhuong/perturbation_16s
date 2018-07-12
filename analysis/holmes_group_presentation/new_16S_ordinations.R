
#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## This is the analog of the depths_plots.R script, but summing over individual
## genes up to the function level. This means we don't have to subset to
## individual species (the data gets smaller), and the resulting names are more
## biologically meaningful.
##
## author: sankaran.kris@gmail.com
## date: 12/31/2017

###############################################################################
## Setup and libraries
###############################################################################
library("tidyverse")
library("phyloseq")
library("PMA")


load("./processed_physeq.rda")

taxtab <- data.frame(tax_table(ps)) %>%
  mutate(Seq_ID = rownames((.)))

# pca.ihs <- prcomp(t(asinh(norm_counts)), scale = FALSE)
# save(list = c("pca.ihs"), file = "./ordinate_16S.rda")

load("ordinate_16S.rda")

sparse_pca <- PMA::SPC(
  scale(t(asinh(norm_counts)), center = TRUE, scale = FALSE), 
  v = pca.ihs$rotation[, 1],
  K = 5, sumabsv = 30)
save(list = c("pca.ihs", "sparse_pca"), file = "./ordinate_16S.rda")
