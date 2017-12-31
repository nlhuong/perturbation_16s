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
#source("annotation.R")

## Define the parser
parser <- arg_parser("Plot MIDAS output")
parser <- add_argument(parser, "--subdir", help = "The subdirectory of data/ containing all the processed data", default = "metagenomic")
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
function_depths <- depths[1:200, ] %>%
  left_join(annotation) %>%
  group_by(ontology, function_id) %>%
  summarise_at(vars(starts_with("M")), sum, na.rm = TRUE)

## can join in interpretable names too
interpretation <- function_interpretation(function_depths$function_id)
function_depths <- function_depths %>%
  left_join(interpretation) %>%
  select(ontology, function_id, interpretation, starts_with("M"))

mfunc <- function_depths %>%
  gather(Meas_ID, value, starts_with("M")) %>%
  left_join(meas %>% select(ends_with("ID"))) %>%
  left_join(samp %>% select(Samp_ID, Subject, ends_with("Interval")))
