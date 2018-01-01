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
library("feather")
library("readxl")
library("argparser")
library("pheatmap")
#source("annotation.R")

## Define the parser
parser <- arg_parser("Plot MIDAS output")
parser <- add_argument(parser, "--subdir", help = "The subdirectory of data/ containing all the processed data", default = "metagenomic")
parser <- add_argument(parser, "--ont", help = "Which ontology to use for annotation. Either 'ec', 'figfam', 'go', or 'path'", default = "go")
argv <- parse_args(parser)

###############################################################################
## Read-in and annotate gene depths with function IDs
###############################################################################
merged_dir <- file.path("..", "data", argv$subdir, "merged")
depths <- read_feather(file.path(merged_dir, "depths.feather"))
meas <- read_xlsx("../data/Mapping_Files_7bDec2017.xlsx", "Meas", skip = 1) %>%
  rename(Samp_ID = SampID)
samp <- read_xlsx("../data/Mapping_Files_7bDec2017.xlsx", "Samp", skip = 1) %>%
  mutate(
    Diet_Interval = factor(Diet_Interval, c("PreDiet", "MidDiet", "PostDiet")),
    CC_Interval = factor(CC_Interval, c("PreCC", "MidCC", "PostCC")),
    Abx_Interval = factor(Abx_Interval, c("PreAbx", "MidAbx", "PostAbx"))
  )

## sum over functions
annotation <- function_annotation(unique(depths$species))
f_depths <- depths[1:10000, ] %>%
  left_join(annotation) %>%
  group_by(ontology, function_id) %>%
  summarise_at(
    vars(starts_with("M")),
    function(x) {
      sum(asinh(x), na.rm = TRUE)
    }
  )

## can join in interpretable names too
interpretation <- function_interpretation(f_depths$function_id)
f_depths <- f_depths %>%
  left_join(interpretation) %>%
  select(ontology, function_id, interpretation, starts_with("M")) %>%
  filter(ontology == argv$ont) %>%
  ungroup()

f_mat <- f_depths %>%
  select_if(is.numeric) %>%
  as.matrix()
rownames(f_mat) <- f_depths$function_id

###############################################################################
## Make a heatmap of these summed function depths
###############################################################################
## remove measurements that are always 0
keep_ix <- colSums(f_mat) > 0
f_mat <- f_mat[, keep_ix]
f_depths <- f_depths %>%
  select_at(!starts_with("M"), colnames(f_mat))

## order functions and measurements by hierarchical clustering
hm <- pheatmap(f_mat, silent = TRUE)
mfunc <- f_depths %>%
  gather(Meas_ID, value, starts_with("M")) %>%
  left_join(meas %>% select(ends_with("ID"))) %>%
  left_join(samp %>% select(Samp_ID, Subject, ends_with("Interval"))) %>%
  mutate(
    Meas_ID = factor(Meas_ID, levels = colnames(f_mat)[hm$tree_col$order]),
    function_id = factor(function_id, levels = rownames(f_mat)[hm$tree_row$order])
  )

ggplot(mfunc) +
  geom_tile(
    aes(x = function_id, y = Meas_ID, alpha = value)
  ) +
  facet_grid(
    Subject ~ function_id, scale = "free", space = "free"
  )
