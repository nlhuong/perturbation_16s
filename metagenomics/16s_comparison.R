#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## The purpose of this script is to compare the abundance of species identified
## through metagenomic vs. 16s profiling. As is, this is a pretty coarse view,
## looking only at species assignments as estimated through MIDAS. A more
## careful analysis would consider the acutal counts of 16s sequences alone.
##
## author: sankaran.kris@gmail.com
## date: 12/11/2017

###############################################################################
## Libraries and setup
###############################################################################
library("tidyverse")
library("phyloseq")

ps <- readRDS("../data/16s/perturb_physeq_8Dec.rds")
metag <- read_tsv("../data/metagenomic/merged/count_reads.tsv") %>%
  separate(species_id, c("genus", "species", "strain_id"), remove = FALSE)

###############################################################################
## Compare genus counts
###############################################################################
tax_table(ps) %>%
  head()

melt_metag <- metag %>%
  gather(sample, count, starts_with("M"))

melt_metag %>%
  group_by(sample, genus) %>%
  summarise(total = sum(count)) %>%
  arrange(sample, desc(total))
