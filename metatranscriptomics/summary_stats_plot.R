#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## The purpose of this script is to generate plots and statistics for
## the the number of reads retained at each step of metatranscriptomics
## data processing.
##
## author: nlhuong90@gmail.com
## date: 03/22/2018

###############################################################################
## Libraries and setup
###############################################################################
library("tidyverse")
library("readxl")
library("phyloseq")
library("ggrepel")

scale_colour_discrete <- function(...)
  scale_colour_brewer(..., palette="Set2")
scale_fill_discrete <- function(...)
  scale_fill_brewer(..., palette="Set2")

theme_set(theme_bw())
theme_updates <- theme(
  panel.border = element_rect(size = 0.5),
  panel.grid = element_blank(),
  axis.ticks = element_blank(),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 6),
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 8),
  strip.background = element_blank(),
  strip.text = element_text(size = 8),
  legend.key = element_blank()
)

###############################################################################
## Data Loading
###############################################################################

# loading sample info
meas <- read_xlsx("../data/Mapping_Files_22Jan2018.xlsx", "Meas", skip = 1) %>%
  rename(Samp_ID = SampID)
interv_levs <- c("NoInterv", "PreDiet", "MidDiet", "PostDiet", "PreCC", "MidCC",
                 "PostCC", "PreAbx", "MidAbx", "PostAbx")
samp <- read_xlsx("../data/Mapping_Files_22Jan2018.xlsx", "Samp", skip = 1) %>%
  mutate(
    Diet_Interval = ifelse(Diet_Interval == "NA", "NoInterv", Diet_Interval),
    CC_Interval = ifelse(CC_Interval == "NA", "NoInterv", CC_Interval),
    Abx_Interval = ifelse(Abx_Interval == "NA", "NoInterv", Abx_Interval),
    Diet_Interval = factor(Diet_Interval, interv_levs),
    CC_Interval = factor(CC_Interval, interv_levs),
    Abx_Interval = factor(Abx_Interval, interv_levs)
  ) %>%
  filter(Samp_Type != "ExtrCont")

# coverage stats
stats <- read_csv("../data/metatranscriptomic/DBU_stats.csv")

stats <- stats %>%
  mutate(
    Subject = strsplit(Meas_ID, "_")[[1]][2], 
    Meas_ID = strsplit(Meas_ID, "_")[[1]][1]) %>%
  select(Meas_ID, Subject, input_fwd, input_rev, 
         trimmed, qual_fltr, unique,
         vector_bwa, vector_blat, human_bwa, human_blat, mRNA,
         # kaiju
         dmnd_refseq, dmnd_seed,
         contigs, bwa_mcds_contigs_ann, dmnd_nr_contigs,
         unassembled, bwa_mcds_unassembled_ann, dmnd_nr_unassembled)

preproc_cols <- c("input_fwd", "input_rev", "trimmed", "qual_fltr", "unique",
                  "vector_bwa", "vector_blat", "human_bwa", "human_blat", "mRNA",
                  "contigs", "unassembled")
melt_preproc <- stats[, c("Meas_ID", preproc_cols)] %>%
  mutate(Meas_ID)
  gather(step, count, input_fwd:unassembled) %>%
  mutate(step = factor(step, levels = preproc_cols))


ggplot(melt_preproc) +
  geom_point(aes(x = step, y=count, color=))

