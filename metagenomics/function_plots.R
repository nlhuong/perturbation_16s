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
  scale_fill_brewer(..., palette = "Set1")
scale_color_interval <- function(...)
  scale_color_brewer(..., palette = "Set1")
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
  mfunc_concat <- mfunc %>%
    gather(perturb_type, perturb, ends_with("Interval")) %>%
    unite(func_pert, function_id, perturb_type, remove = FALSE)
  pert_types <- unique(mfunc_concat$perturb_type)
  func_pert_lev <-  paste(
    rep(levels(mfunc$function_id), times = length(pert_types)),
    rep(pert_types, each = nlevels(mfunc$function_id)),
    sep = "_"
  )
  mfunc_concat$func_pert <- factor(
    mfunc_concat$func_pert,
    levels = func_pert_lev
  )

  ggplot(mfunc_concat) +
    geom_tile(
      aes(x = func_pert, y = Meas_ID, alpha = sqrt(value), fill = perturb)
    ) +
    facet_grid(Subject ~ ., scale = "free", space = "free") +
    scale_alpha(range = c(0.01, 1)) +
    scale_fill_interval() +
    theme(
      axis.text = element_blank(),
      legend.position = "bottom"
    )
}

###############################################################################
## Read-in and annotate gene depths with function IDs
###############################################################################
merged_dir <- file.path("..", "data", argv$subdir, "merged")
depths <- read_feather(file.path(merged_dir, "depths.feather"))
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

meas <- read_xlsx("../data/Mapping_Files_7bDec2017.xlsx", "Meas", skip = 1) %>%
  rename(Samp_ID = SampID) %>%
  left_join(samp)
bad_samples <- c("M3728", "M3673", "M3695", "M3204", "M3064", "M3109", "M3188",
                 "M3654")

## sum over function IDs (after asinh transforming)
## annotation <- function_annotation(unique(depths$species))
## write.csv(annotation, "annotation.csv", row.names = FALSE)
annotation <- read.csv("annotation.csv")
f_depths <- depths %>%
  left_join(annotation) %>%
  group_by(ontology, function_id) %>%
  summarise_at(
    vars(starts_with("M")),
    function(x) {
      sum(asinh(x), na.rm = TRUE)
    }
  )
f_depths <- f_depths[, !(colnames(f_depths) %in% bad_samples)]

## can join in interpretable names too
## interpretation <- function_interpretation(unique(annotation$function_id))
## write.csv(interpretation, "interpretation.csv", row.names = FALSE)
interpretation <- read.csv("interpretation.csv")
f_depths <- f_depths %>%
  left_join(interpretation) %>%
  select(ontology, function_id, interpretation, starts_with("M")) %>%
  filter(ontology == argv$ont) %>%
  ungroup()

f_mat <- f_depths %>%
  ungroup() %>%
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
## mfunc <- melt_fmat(f_depths, f_mat, meas, samp)
## plot_mfunc(mfunc)
## ggsave("go_heatmap.png", width = 13.4, height = 5.4)

###############################################################################
## Same code as above, but after variance filtering the GO IDs
###############################################################################
keep_ix <- apply(log(1 + f_mat), 1, var) > 0.1
f_mat <- f_mat[keep_ix, ]
f_depths <- f_depths[keep_ix, ]

mfunc <- melt_fmat(f_depths, f_mat, meas, samp)
plot_mfunc(mfunc) +
  labs(
    x = "Function",
    y = "Sample",
    col = "State",
    alpha = "sqrt(coverage)"
  )
ggsave("go_heatmap_var_filter.png", dpi = 500, width = 8.64, height = 4.04)

###############################################################################
## Sparse pca
###############################################################################
pc_func <- SPC(scale(f_mat), K = 5, sumabsv = 18)
scores <- data.frame(
  "function_id" = rownames(f_mat),
  pc_func$u %*% diag(pc_func$d)
) %>%
  left_join(f_depths %>% select(-starts_with("M")))

loadings <- data.frame(
  pc_func$v,
  "Meas_ID" = colnames(f_mat)
) %>%
  left_join(meas %>% select(ends_with("ID"))) %>%
  left_join(samp %>% select(Samp_ID, Subject, ends_with("Interval")))

d1 <- 100 * pc_func$d[1] / sum(pc_func$d)
d2 <- 100 * pc_func$d[2] / sum(pc_func$d)
ggplot(scores) +
  geom_hline(yintercept = 0, col = "#e6e6e6") +
  geom_vline(xintercept = 0, col = "#e6e6e6") +
  geom_abline(slope = -1, col = "#e6e6e6") +
  geom_point(
    aes(x = X1, y = X2),
    size = 1,
    alpha = 0.4
  ) +
  labs(
    x = sprintf("Axis 1 [%s%%]", round(d1, 3)),
    y = sprintf("Axis 2 [%s%%]", round(d2, 3))
  ) +
  coord_fixed(sqrt(d2 / d1)) +
  geom_text_repel(
    data = scores %>%
      filter(abs(X1 - X2) > 80),
    aes(x = X1, y = X2, label = interpretation),
    size = 2,
    force = 0.1
  ) +
  theme(legend.position = "bottom")
ggsave("function_scores.png", dpi = 500, width = 6.17, height = 1.81)

ggplot(loadings) +
  geom_point(
    aes(x = X1, y = X2, col = Diet_Interval),
    size = 1,
    pch = 15,
    position = position_nudge(y = -0.006),
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
    position = position_nudge(y = 0.006),
    alpha = 0.6
  ) +
  labs(
    x = sprintf("Axis 1 [%s%%]", round(d1, 3)),
    y = sprintf("Axis 2 [%s%%]", round(d2, 3))
  ) +
  scale_color_interval() +
  guides(col = guide_legend(override.aes = list(size = 6))) +
  facet_wrap(~ Subject, scales = "free") +
  theme(legend.position = "bottom")
ggsave("function_loadings.png", dpi = 500, width = 6.17, height = 3.15)
