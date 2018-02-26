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
library("argparser")
library("readxl")
library("feather")
library("ggrepel")
library("pheatmap")
library("forcats")
library("PMA")
library("tidyverse")

## Define the parser
parser <- arg_parser("Plot MIDAS output")
parser <- add_argument(parser, "--subdir", help = "The subdirectory of data/ containing all the processed data", default = "metagenomic")
parser <- add_argument(parser, "--k", help = "k in k-over-a filter for coverage", default = 0.05)
parser <- add_argument(parser, "--a", help = "a in k-over-a filter for coverage", default = 0.07)
argv <- parse_args(parser)

## custom plot defaults
scale_fill_interval <- function(...)
  scale_fill_brewer(..., palette = "Set1")
scale_color_interval <- function(...)
  scale_color_brewer(..., palette = "Set1")
scale_colour_discrete <- function(...)
  scale_colour_brewer(..., palette="Set2")

theme_set(theme_bw())
theme_update(
  panel.border = element_rect(size = 0.5),
  panel.background = element_rect(fill = "#f7f7f7"),
  panel.grid = element_blank(),
  axis.ticks = element_blank(),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 6),
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 8),
  strip.background = element_blank(),
  strip.text = element_text(size = 8),
  legend.key = element_blank(),
  legend.background = element_rect(fill = "#dadada")
)

## Functions used throughout
plot_cov <- function(mcoverage) {
  mcoverage_concat <- mcoverage %>%
    gather(perturb_type, perturb, ends_with("Interval")) %>%
    unite(species_pert, species_id, perturb_type, remove = FALSE)
  pert_types <- unique(mcoverage_concat$perturb_type)
  species_pert_lev <-  paste(
    rep(levels(mcoverage$species_id), times = length(pert_types)),
    rep(pert_types, each = nlevels(mcoverage$species_id)),
    sep = "_"
  )
  mcoverage_concat$species_pert <- factor(
    mcoverage_concat$species_pert,
    levels = species_pert_lev
  )

  ggplot(mcoverage_concat) +
    geom_tile(
      aes(x = species_pert, y = Meas_ID, alpha = coverage, fill = perturb)
    ) +
    facet_grid(
      Subject ~ genus_grouped,
      scales = "free",
      space = "free"
    ) +
    scale_fill_interval() +
    scale_alpha(range = c(0.01, 1)) +
    theme(
      axis.text = element_blank(),
      axis.title.y = element_blank(),
      strip.text.x = element_text(angle = 90, hjust = 0),
      strip.text.y = element_text(angle = 0),
      legend.position = "bottom"
    )
}

plot_loadings <- function(loadings, d1 = 1, d2 = 1) {
  ggplot(loadings) +
    geom_point(
      aes(x = X1, y = X2, col = Diet_Interval),
      size = 1,
      pch = 15,
      position = position_nudge(y = -0.0015),
      alpha = 0.6
    ) +
    geom_point(
      aes(x = X1, y = X2, col = CC_Interval),
      size = 1,
      pch = 15,
      alpha = 0.6
    ) +
    geom_point(
      aes(x = X1, y = X2, col = Abx_Interval),
      size = 1,
      pch = 15,
      position = position_nudge(y = 0.0015),
      alpha = 0.6
    ) +
    labs(
      x = sprintf("Axis 1 [%s%%]", round(d1, 3)),
      y = sprintf("Axis 2 [%s%%]", round(d2, 3))
    ) +
    scale_color_interval() +
    guides(col = guide_legend(override.aes = list(size = 6))) +
    facet_wrap(~ Subject, scales = "free")
}

## Read in data
merged_dir <- file.path("..", "data", argv$subdir, "merged")
coverage <- read_tsv(file.path(merged_dir, "coverage.tsv"))
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
## Prepare data for visualizing species COG data
###############################################################################
keep_ix <- rowMeans(coverage[, -1] > argv$k) > argv$a
coverage <- coverage[keep_ix, ]
coverage[, -1] <- asinh(coverage[, -1])

# some samples seem missing...
meas %>%
  filter(Meas_ID == "M3654")
samp %>%
  filter(Samp_ID == "S1639")
bad_samples <- c("M3728", "M3673", "M3695", "M3204", "M3064", "M3109", "M3188",
                 "M3654")
coverage <- coverage[, !(colnames(coverage) %in% bad_samples)]

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
  gather(Meas_ID, coverage,  starts_with("M")) %>%
  left_join(meas %>% select(Meas_ID, Samp_ID)) %>%
  left_join(samp %>% select(Samp_ID, Subject, Samp_Date, ends_with("Interval")))

rownames(coverage) <- coverage$species_id
cov_mat <- coverage %>%
  select(starts_with("M"))
sp_order <- pheatmap(cov_mat, silent = TRUE, cluster_cols = FALSE)$tree_row$order
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
  dplyr::summarise(sum = sum(coverage)) %>%
  arrange(desc(sum)) %>%
  .[["genus_grouped"]]
mcoverage <- mcoverage %>%
  mutate(
    Meas_ID = factor(Meas_ID, meas_levels),
    genus_grouped = factor(genus_grouped, genus_levels)
  )

## pca biplot data
pc_cov <- SPC(scale(cov_mat), K = 5, sumabsv = 15)
scores <- data.frame(
  "species_id" = rownames(cov_mat),
  pc_cov$u %*% diag(pc_cov$d)
) %>%
  left_join(coverage %>% select(-starts_with("M")))

loadings <- data.frame(
  pc_cov$v,
  "Meas_ID" = colnames(cov_mat)
) %>%
  left_join(meas %>% select(ends_with("ID"))) %>%
  left_join(samp %>% select(Samp_ID, Subject, ends_with("Interval")))

###############################################################################
## Produce visualizations for the species COG data
###############################################################################
plot_cov(mcoverage)
ggsave("coverage_heatmap.png", dpi = 500)

## basic pca biplot
d1 <- 100 * pc_cov$d[1] / sum(pc_cov$d)
d2 <- 100 * pc_cov$d[2] / sum(pc_cov$d)
ggplot(scores) +
  geom_hline(yintercept = 0, col = "#e6e6e6") +
  geom_vline(xintercept = 0, col = "#e6e6e6") +
  geom_point(
    aes(x = X1, y = X2, col = genus_grouped),
    size = 1
  ) +
  labs(
    x = sprintf("Axis 1 [%s%%]", round(d1, 3)),
    y = sprintf("Axis 2 [%s%%]", round(d2, 3))
  ) +
  coord_fixed(sqrt(d2 / d1)) +
  geom_text_repel(
    data = scores %>%
      filter(X1 ^ 2 + X2 ^ 2 > 100),
    aes(x = X1, y = X2, label = species, col = genus_grouped)
  ) +
  theme(legend.position = "bottom")
ggsave("coverage_scores.png", height = 5.09, width = 5.09, dpi = 500)

## presumably the two subjects
plot_loadings(loadings, d1, d2) +
  theme(legend.position = "bottom")
ggsave("coverage_loadings.png", width = 7.89, height = 4.04, dpi = 500)
