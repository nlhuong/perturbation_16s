#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## A first pass at visualizing the metagenomics data produced by using MIDAS
## with the default settings.
##
## author: sankaran.kris@gmail.com
## date: 12/20/2017

###############################################################################
## Setup and data
###############################################################################
library("tidyverse")
library("argparser")
library("readxl")
library("feather")
library("pheatmap")
library("ggrepel")
library("forcats")

## Define the parser
parser <- arg_parser("Plot MIDAS output")
parser <- add_argument(parser, "--subdir", help = "The subdirectory of data/ containing all the processed data", default = "metagenomic")
parser <- add_argument(parser, "--k", help = "k in k-over-a filter for depths", default = 0.25)
parser <- add_argument(parser, "--a", help = "a in k-over-a filter for depths", default = 2)
argv <- parse_args(parser)

## custom plot defaults
scale_fill_interval <- function(...)
  scale_fill_manual(..., values = c("#15c7b0", "#c71585", "#15c7b0"), na.value = "#2f4f4f")
scale_color_interval <- function(...)
  scale_color_manual(..., values = c("#15c7b0", "#c71585", "#15c7b0"), na.value = "#2f4f4f")
scale_colour_discrete <- function(...)
  scale_colour_brewer(..., palette = "Set2")

theme_set(theme_bw())
theme_update(
  panel.border = element_rect(size = 0.5),
  panel.grid = element_blank(),
  axis.ticks = element_blank(),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 9),
  axis.text = element_text(size = 9),
  axis.title = element_text(size = 8),
  strip.background = element_blank(),
  strip.text = element_text(size = 9),
  legend.key = element_blank()
)

## Functions used throughout
plot_heatmap <- function(mdepths, fill_type = "Diet_Interval") {
  ggplot(mdepths) +
    geom_tile(
      aes_string(
        x = "gene_id",
        y = "Meas_ID",
        alpha = "depth",
        fill = fill_type
      )
    ) +
    scale_alpha(range = c(0, 1)) +
    scale_fill_interval() +
    facet_grid(Subject ~ species_name, scale = "free", space = "free") +
    theme(
      axis.text = element_blank(),
      strip.text.y = element_text(angle = 0),
      legend.position = "bottom"
    )
}

plot_loadings <- function(loadings, col_type) {
  ggplot(loadings) +
    geom_hline(yintercept = 0, col = "#e6e6e6") +
    geom_vline(xintercept = 0, col = "#e6e6e6") +
    geom_text_repel(
      aes_string(x = "Comp.1", y = "Comp.2", label = "Meas_ID", col = col_type),
      size = 2,
      force = 0.001
    ) +
    guides(col = guide_legend(override.aes = list(size = 6))) +
    scale_color_interval() +
    facet_grid(~ Subject)
}

## Read in data
merged_dir <- file.path("..", "data", argv$subdir, "merged")
depths <- read_feather(file.path(merged_dir, "depths.feather")) %>%
  mutate(genus = str_extract(species, "[^_]+"))
meas <- read_xlsx("../data/Mapping_Files_7bDec2017.xlsx", "Meas", skip = 1) %>%
  rename(Samp_ID = SampID)
interv_levs <- c("NoInterv", "PreDiet", "MidDiet", "PostDiet", "PreCC", "MidCC",
                 "PostCC", "PreAbx", "MidAbx", "PostAbx")
samp <- read_xlsx("../data/Mapping_Files_7bDec2017.xlsx", "Samp", skip = 1) %>%
  mutate(
    Diet_Interval = ifelse(Diet_Interval == "NA", "NoInterv", Diet_Interval),
    CC_Interval = ifelse(CC_Interval == "NA", "NoInterv", CC_Interval),
    Abx_Interval = ifelse(Abx_Interval == "NA", "NoInterv", Abx_Interval),
    Diet_Interval = factor(Diet_Interval, interv_levs),
    CC_Interval = factor(CC_Interval, interv_levs),
    Abx_Interval = factor(Abx_Interval, interv_levs)
  ) %>%
  filter(Samp_Type != "ExtrCont")

###############################################################################
## Prepare data for gene depths plots
###############################################################################
depths <- depths %>%
  filter(genus %in% c("Parabacteroides", "Eubacterium", "Bacteroides"))
depths[is.na(depths)] <- 0
keep_ix <- rowMeans(depths[, -c(1, 2)] > argv$a) > argv$k

depths_df <- depths[keep_ix, ] %>%
  as.data.frame() %>%
  select(-genus) %>%
  separate(
    species,
    c("genus", "species_name", "strain_id"), "_",
    remove = FALSE
  ) %>%
  mutate_at(vars(starts_with("M")), asinh)
rownames(depths_df) <- depths_df$gene_id
rm(depths)

## transformed depths
depths_df[1:10000, ] %>%
  select(starts_with("M")) %>%
  as.matrix() %>%
  hist(breaks = 30, main = "asinh(Depths)", ylim = c(0, 2e5))

## extracting order for plot
meas_levels <- meas %>%
  left_join(samp) %>%
  select(Meas_ID, Subject, Samp_Date) %>%
  unique() %>%
  arrange(Subject, desc(Samp_Date)) %>%
  .[["Meas_ID"]]
depths_mat <- depths_df %>%
  select(starts_with("M")) %>%
  as.matrix()
col_sums <- colSums(depths_mat)
bad_ix <- col_sums < 8e3 | col_sums > 5.5e4
depths_mat <- depths_mat[, !bad_ix]

hm <- pheatmap(depths_mat, silent = TRUE, cluster_cols = FALSE)

mdepths <- depths_df %>%
  gather(Meas_ID, depth, starts_with("M")) %>%
  left_join(meas %>% select(ends_with("ID"))) %>%
  left_join(samp %>% select(Samp_ID, Subject, ends_with("Interval"))) %>%
  mutate(
    Meas_ID = factor(Meas_ID, levels = meas_levels),
    gene_id = factor(gene_id, levels = rownames(depths_mat)[hm$tree_row$order])
  )

## data for PCA biplot
depths_mat <- depths_mat[, apply(depths_mat, 2, var) > 0]
pc_depths <- princomp(scale(depths_mat, scale = FALSE))
scores <- pc_depths$scores[, 1:10] %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  left_join(depths_df %>% select(gene_id, species))

loadings <- pc_depths$loadings[, 1:5] %>%
  as.data.frame() %>%
  rownames_to_column("Meas_ID") %>%
  left_join(meas %>% select(ends_with("ID"))) %>%
  left_join(samp %>% select(Samp_ID, Subject, ends_with("Interval")))

###############################################################################
## Make plots for the analysis
###############################################################################
plot_heatmap(mdepths, "Diet_Interval")
plot_heatmap(mdepths, "CC_Interval")
plot_heatmap(mdepths, "Abx_Interval")

ggplot(scores) +
  geom_hline(yintercept = 0, col = "#e6e6e6") +
  geom_vline(xintercept = 0, col = "#e6e6e6") +
  geom_point(
    aes(x = Comp.1, y = Comp.2, size = Comp.3, col = species),
    alpha = 0.2
  ) +
  geom_text_repel(
    data = scores %>%
      filter(Comp.1 ^ 2 + Comp.2 ^ 2 > 3500),
    aes(x = Comp.1, y = Comp.2, label = gene_id),
    force = 0.2,
    size = 2.5
  ) +
  guides(
    col = guide_legend(override.aes = list(alpha = 1)),
    size = guide_legend(override.aes = list(alpha = 1))
  ) +
  scale_size(range = c(0.001, 2))

plot_loadings(loadings, "Diet_Interval")
plot_loadings(loadings, "CC_Interval")
plot_loadings(loadings, "Abx_Interval")
