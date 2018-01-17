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
source("annotation.R")
Sys.setenv("MIDAS_DB" = "/scratch/users/kriss1/applications/MIDAS/database/midas_db_v1.2")

scale_fill_interval <- function(...)
  scale_fill_manual(..., values = c("#15c7b0", "#c71585", "#15c7b0"), na.value = "#2f4f4f")
scale_colour_discrete <- function(...)
  scale_colour_brewer(..., palette="Set2")
scale_fill_discrete <- function(...)
  scale_fill_brewer(..., palette="Set2")

theme_set(theme_bw())
theme_update(
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

## Define the parser
parser <- arg_parser("Plot MIDAS output")
parser <- add_argument(parser, "--subdir", help = "The subdirectory of data/ containing all the processed data", default = "metagenomic")
parser <- add_argument(parser, "--ont", help = "Which ontology to use for annotation. Either 'ec', 'figfam', 'go', or 'path'", default = "go")
argv <- parse_args(parser)

###############################################################################
## Helper functions
###############################################################################
#' Melt depths, with heatmap factor level ordering
#'
#' This melts the function abundance matrix. We could just use pheatmap, but
#' then we wouldn't be able to facet by subject, say.
melt_fmat <- function(f_depths, f_mat, meas, samp) {
  hm <- pheatmap(f_mat, silent = TRUE)
  meas_levels <- meas %>%
    select(Meas_ID, Subject, Samp_Date) %>%
    unique() %>%
    arrange(Subject, desc(Samp_Date)) %>%
    .[["Meas_ID"]]

  f_depths %>%
    gather(Meas_ID, value, starts_with("M")) %>%
    left_join(meas %>% select(ends_with("ID"))) %>%
    left_join(samp %>% select(Samp_ID, Subject, ends_with("Interval"))) %>%
    mutate(
      Meas_ID = factor(Meas_ID, levels = meas_levels),
      function_id = factor(function_id, levels = rownames(f_mat)[hm$tree_row$order])
    )
}

#' Plot the Melted Function IDs
plot_mfunc <- function(mfunc) {
  ggplot(mfunc) +
    geom_tile(
      aes(x = function_id, y = Meas_ID, alpha = log(1 + value), fill = Diet_Interval)
    ) +
    facet_grid(Subject ~ ., scale = "free", space = "free") +
    scale_alpha(range = c(0, 1)) +
    scale_fill_interval() +
    theme(
      axis.text.y = element_text(size = 4),
      axis.text.x = element_text(size = 3, angle = -90)
    )
}

###############################################################################
## Read-in and annotate gene depths with function IDs
###############################################################################
merged_dir <- file.path("..", "data", argv$subdir, "merged")
depths <- read_feather(file.path(merged_dir, "depths.feather"))
samp <- read_xlsx("../data/Mapping_Files_7bDec2017.xlsx", "Samp", skip = 1) %>%
  mutate(
    Diet_Interval = factor(Diet_Interval, c("PreDiet", "MidDiet", "PostDiet")),
    CC_Interval = factor(CC_Interval, c("PreCC", "MidCC", "PostCC")),
    Abx_Interval = factor(Abx_Interval, c("PreAbx", "MidAbx", "PostAbx"))
  )
meas <- read_xlsx("../data/Mapping_Files_7bDec2017.xlsx", "Meas", skip = 1) %>%
  rename(Samp_ID = SampID) %>%
  left_join(samp)

## sum over function IDs (after asinh transforming)
annotation <- function_annotation(unique(depths$species))
f_depths <- depths %>%
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
## remove measurements that are always 0 (why does this happen?)
keep_ix <- colSums(f_mat) > 0
f_mat <- f_mat[, keep_ix]
f_depths <- f_depths %>%
  select_at(vars(-starts_with("M"), colnames(f_mat)))

## order functions and measurements by hierarchical clustering
mfunc <- melt_fmat(f_depths, f_mat, meas, samp)
plot_mfunc(mfunc)
ggsave("go_heatmap.png", width = 13.4, height = 5.4)

###############################################################################
## Same code as above, but after variance filtering the GO IDs
###############################################################################
keep_ix <- apply(log(1 + f_mat), 1, var) > 0.3
f_mat <- f_mat[keep_ix, ]
f_depths <- f_depths[keep_ix, ]

mfunc <- melt_fmat(f_depths, f_mat, meas, samp)
plot_mfunc(mfunc)
ggsave("go_heatmap_var_filter.png", width = 13.4, height = 5.4)
