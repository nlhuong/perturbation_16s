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
library("argparser")
library("readxl")
library("feather")
library("pheatmap")
library("ggrepel")
library("forcats")

## custom plot defaults
scale_colour_discrete <- function(...)
  scale_colour_brewer(..., palette="Set2")
scale_fill_discrete <- function(...)
  scale_fill_brewer(..., palette="Set2")
scale_fill_interval <- function(...)
  scale_fill_manual(..., values = c("#15c7b0", "#c71585", "#15c7b0", "red"), na.value = "#2f4f4f")

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
plot_cov <- function(hm, cov_aes, fill_type) {
  cov_aes$fill <- fill_type
  hm +
    geom_tile(cov_aes)
}

## Define the parser
parser <- arg_parser("Plot MIDAS output")
parser <- add_argument(parser, "--subdir", help = "The subdirectory of data/ containing all the processed data", default = "metagenomic")
parser <- add_argument(parser, "--k_cov", help = "k in k-over-a filter for coverage", default = 0.05)
parser <- add_argument(parser, "--a_cov", help = "a in k-over-a filter for coverage", default = 0)
parser <- add_argument(parser, "--k_depth", help = "k in k-over-a filter for depths", default = 0.1)
parser <- add_argument(parser, "--a_depth", help = "a in k-over-a filter for depths", default = 0.0)
argv <- parse_args(parser)

## Read in data
merged_dir <- file.path("..", "data", argv$subdir, "merged")
coverage <- read_tsv(file.path(merged_dir, "coverage.tsv"))
depths <- read_feather(file.path(merged_dir, "depths.feather"))
copy_num <- read_feather(file.path(merged_dir, "copy_num.feather"))
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
keep_ix <- rowMeans(coverage[, -1] > argv$k_cov) > argv$a_cov
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
    genus_grouped = fct_lump(genus, prop = 0.04)
  )

mcoverage <- coverage %>%
  gather(Meas_ID, coverage, starts_with("M")) %>%
  left_join(meas %>% select(Meas_ID, Samp_ID)) %>%
  left_join(samp %>% select(Samp_ID, Subject, Samp_Date, ends_with("Interval")))

rownames(coverage) <- coverage$species_id
cov_mat <- coverage %>%
  select(starts_with("M"))
sp_order <- hclust(dist(cov_mat))$order
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
pc_cov <- princomp(scale(t(cov_mat)))
scores <- pc_cov$scores[, 1:5] %>%
  as.data.frame() %>%
  rownames_to_column("species_id") %>%
  left_join(coverage %>% select(-starts_with("M")))

loadings <- pc_cov$loadings[, 1:5] %>%
  as.data.frame() %>%
  rownames_to_column("sample")


###############################################################################
## Produce visualizations for the species COG data
###############################################################################
cov_aes <- aes(
  x = species_id,
  y = Meas_ID,
  alpha = coverage,
  fill = Diet_Interval
)
hm <- ggplot(mcoverage) +
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

plot_cov(hm, cov_aes, substitute(CC_Interval))
plot_cov(hm, cov_aes, substitute(Diet_Interval))
plot_cov(hm, cov_aes, substitute(Abx_Interval))

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
ggplot(loadings) +
  geom_text_repel(
    aes(x = Comp.1, y = Comp.2, label = sample),
    size = 2,
    force = 0.001
  )

###############################################################################
## Gene depths
###############################################################################
depths[is.na(depths)] <- 0
keep_ix <- rowMeans(depths[, -c(1, 2)] > argv$a_depth) > argv$k_depth

depths_df <- depths[keep_ix, ] %>%
  as.data.frame() %>%
  separate(
    species,
    c("genus", "species_name", "strain_id"), "_",
    remove = FALSE
  ) %>%
  mutate(
    genus = as.factor(genus),
    genus_grouped = fct_lump(genus, prop = 0.04)
  ) %>%
  mutate_at(
    vars(starts_with("M")),
    .funs = function(x) {
      x[is.na(x)] <- 0
      asinh(x)
    })
rm(depths)

## transformed depths
depths_df[1:1000, ] %>%
  select(starts_with("M")) %>%
  as.matrix() %>%
  hist(breaks = 20)

pheatmap(
  depths_df %>%
  filter(genus == "Ruminococcus") %>%
  select(starts_with("M")) %>%
  as.matrix() %>%
  t(),
  show_rownames = FALSE
)

## can also make a PCA biplot
depths_mat <- depths_df %>%
  filter(genus == "Ruminococcus") %>%
  select(starts_with("M")) %>%
  as.matrix()
rownames(depths_mat) <- depths_df %>%
  filter(genus == "Ruminococcus") %>%
  .[["gene_id"]]

pc_depths <- princomp(depths_mat)
scores <- pc_depths$scores %>%
  scale() %>%
  as.data.frame() %>%
  rownames_to_column("gene_id")

ggplot(scores) +
  geom_hline(yintercept = 0, col = "#e6e6e6") +
  geom_vline(xintercept = 0, col = "#e6e6e6") +
  geom_point(
    aes(x = Comp.1, y = Comp.2, size = Comp.3),
    alpha = 0.4
  ) +
  geom_text_repel(
    data = scores %>%
      filter(Comp.1 ^ 2 + Comp.2 ^ 2 > 18),
    aes(x = Comp.1, y = Comp.2, label = gene_id),
    force = 0.002,
    size = 2.5
  ) +
    scale_size(range = c(0.001, 2)) +
    ylim(-4, 8)
