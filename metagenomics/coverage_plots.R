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
library("ggrepel")
library("forcats")

## custom plot defaults
scale_fill_interval <- function(...)
  scale_fill_manual(..., values = c("#15c7b0", "#c71585", "#15c7b0"), na.value = "#2f4f4f")
scale_color_interval <- function(...)
  scale_color_manual(..., values = c("#15c7b0", "#c71585", "#15c7b0"), na.value = "#2f4f4f")
scale_colour_discrete <- function(...)
  scale_colour_brewer(..., palette="Set2")

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
plot_cov <- function(mcoverage, fill_type) {
  ggplot(mcoverage) +
    geom_tile(
      aes_string(x = "species_id", y = "Meas_ID", alpha = "coverage", fill = fill_type)
    ) +
    facet_grid(
      Subject ~ genus_grouped,
      scales = "free",
      space = "free"
    ) +
    scale_fill_interval() +
    scale_alpha(range = c(0, 1)) +
    theme(
      axis.text = element_blank(),
      axis.title.y = element_blank(),
      strip.text.x = element_text(angle = 90, hjust = 0),
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

## Define the parser
parser <- arg_parser("Plot MIDAS output")
parser <- add_argument(parser, "--subdir", help = "The subdirectory of data/ containing all the processed data", default = "metagenomic")
parser <- add_argument(parser, "--k", help = "k in k-over-a filter for coverage", default = 0.05)
parser <- add_argument(parser, "--a", help = "a in k-over-a filter for coverage", default = 0)
argv <- parse_args(parser)

## Read in data
merged_dir <- file.path("..", "data", argv$subdir, "merged")
coverage <- read_tsv(file.path(merged_dir, "coverage.tsv"))
meas <- read_xlsx("../data/Mapping_Files_7bDec2017.xlsx", "Meas", skip = 1) %>%
  rename(Samp_ID = SampID)
samp <- read_xlsx("../data/Mapping_Files_7bDec2017.xlsx", "Samp", skip = 1) %>%
  mutate(
    Diet_Interval = factor(Diet_Interval, c("PreDiet", "MidDiet", "PostDiet")),
    CC_Interval = factor(CC_Interval, c("PreCC", "MidCC", "PostCC")),
    Abx_Interval = factor(Abx_Interval, c("PreAbx", "MidAbx", "PostAbx"))
  )

###############################################################################
## Prepare data for visualizing species COG data
###############################################################################
keep_ix <- rowMeans(coverage[, -1] > argv$k) > argv$a
coverage <- coverage[keep_ix, ]
coverage[, -1] <- asinh(coverage[, -1])

coverage <- as.data.frame(coverage) %>%
  separate(
    species_id,
    c("genus", "species", "strain_id"), "_",
    remove = FALSE
  ) %>%
  mutate(
    genus = as.factor(genus),
    genus_grouped = fct_lump(genus, n = 7)
  )

mcoverage <- coverage %>%
  gather(Meas_ID, coverage, starts_with("M")) %>%
  left_join(meas %>% select(Meas_ID, Samp_ID)) %>%
  left_join(samp %>% select(Samp_ID, Subject, Samp_Date, ends_with("Interval")))

rownames(coverage) <- coverage$species_id
cov_mat <- coverage %>%
  select(starts_with("M"))
sp_order <- pheatmap(cov_mat, silent = TRUE)$tree_row$order
mcoverage$species_id <- factor(
  mcoverage$species_id,
  levels = rownames(cov_mat)[sp_order]
)

## sort measurement and genus levels
meas_levels <- mcoverage %>%
  select(Meas_ID, Subject, Samp_Date) %>%
  unique() %>%
  arrange(Subject, desc(Samp_Date)) %>%
  .[["Meas_ID"]]
genus_levels <- mcoverage %>%
  group_by(genus_grouped) %>%
  summarise(sum = sum(coverage)) %>%
  arrange(desc(sum)) %>%
  .[["genus_grouped"]]
mcoverage <- mcoverage %>%
  mutate(
    Meas_ID = factor(Meas_ID, meas_levels),
    genus_grouped = factor(genus_grouped, genus_levels)
  )

## pca biplot data
pc_cov <- princomp(scale(cov_mat))
scores <- pc_cov$scores[, 1:5] %>%
  as.data.frame() %>%
  rownames_to_column("species_id") %>%
  left_join(coverage %>% select(-starts_with("M")))

loadings <- pc_cov$loadings[, 1:5] %>%
  as.data.frame() %>%
  rownames_to_column("Meas_ID") %>%
  left_join(meas %>% select(ends_with("ID"))) %>%
  left_join(samp %>% select(Samp_ID, Subject, ends_with("Interval")))


###############################################################################
## Produce visualizations for the species COG data
###############################################################################
plot_cov(mcoverage, "CC_Interval")
plot_cov(mcoverage, "Diet_Interval")
plot_cov(mcoverage, "Abx_Interval")

## basic pca biplot
ggplot(scores) +
  geom_hline(yintercept = 0, col = "#e6e6e6") +
  geom_vline(xintercept = 0, col = "#e6e6e6") +
  geom_point(
    aes(x = Comp.1, y = Comp.2, col = genus_grouped),
    size = 1
  ) +
  geom_text_repel(
    data = scores %>%
      filter(Comp.1 ^ 2 + Comp.2 ^ 2 > 100),
    aes(x = Comp.1, y = Comp.2, label = species, col = genus_grouped)
  )

## presumably the two subjects
plot_loadings(loadings, "CC_Interval")
plot_loadings(loadings, "Diet_Interval")
plot_loadings(loadings, "Abx_Interval")
