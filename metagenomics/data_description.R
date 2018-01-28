#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Some characterization of the metagenomic samples before any profiling has
## been done.
##
## author: sankaran.kris@gmail.com
## date: 01/27/2017

library("readxl")
library("tidyverse")
library("RColorBrewer")

scale_colour_discrete <- function(...)
  scale_colour_brewer(..., palette="Set1")
scale_fill_discrete <- function(...)
  scale_fill_brewer(..., palette="Set1")

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

## number of metagenomic samples (considering forwards vs. reversed as one)
fqs <- list.files("../data/metagenomic/", "*.fq", full.names = TRUE)
fq_ids <- str_extract(fqs, "M[0-9]+")
length(fqs) / 2

## metagenomics version of experiment design
interv_levs <- c("NoInterv", "PreDiet", "MidDiet", "PostDiet", "PreCC", "MidCC",
                 "PostCC", "PreAbx", "MidAbx", "PostAbx")
samp <- read_xlsx("../data/Mapping_Files_7bDec2017.xlsx", "Samp", skip = 1) %>%
  mutate(
    Diet_Interval = ifelse(Diet_Interval == "NA", "NoInterv", Diet_Interval),
    CC_Interval = ifelse(CC_Interval == "NA", "NoInterv", CC_Interval),
    Abx_Interval = ifelse(Abx_Interval == "NA", "NoInterv", Abx_Interval),
    Diet_Interval = factor(Diet_Interval, interv_levs),
    CC_Interval = factor(CC_Interval, interv_levs),
    Abx_Interval = factor(Abx_Interval, interv_levs),
  ) %>%
  filter(Samp_Type != "ExtrCont")

meas <- read_xlsx("../data/Mapping_Files_7bDec2017.xlsx", "Meas", skip = 1) %>%
  rename(Samp_ID = SampID)
mg_samp_ids <- meas %>%
  filter(Meas_Type == "MetaG") %>%
  .[["Samp_ID"]]
no_interv <- samp %>%
  group_by(Subject) %>%
  summarise(
    no_interv = all(Diet_Interval == "NoInterv") &&
    all(CC_Interval == "NoInterv") &&
    all(Abx_Interval == "NoInterv")
  ) %>%
  filter(no_interv) %>%
  .[["Subject"]]
sequenced_ids <- meas %>%
  filter(Meas_ID %in% fq_ids) %>%
  .[["Samp_ID"]]
samp$metag_complete <- "none"
samp$metag_complete[samp$Samp_ID %in% mg_samp_ids] <- "collected"
samp$metag_complete[samp$Samp_ID %in% sequenced_ids] <- "sequenced"

samp <- samp %>%
  mutate(
    metag_complete = factor(metag_complete, levels = c("sequenced", "collected", "none")),
    no_interv = Subject %in% no_interv
  )

exp_design_plot <- function(df) {
  ggplot(df) +
    geom_point(
      aes(x = Samp_Date, y = Subject, col = Diet_Interval),
      pch = 73, size = 0.6, position = position_nudge(y = 0.3)
    ) +
    geom_point(
      aes(x = Samp_Date, y = Subject, col = CC_Interval),
      pch = 73, size = 0.6
    ) +
    geom_point(
      aes(x = Samp_Date, y = Subject, col = Abx_Interval),
      pch = 73, size = 0.6, position = position_nudge(y = -0.3)
    ) +
    guides(
      col = guide_legend(override.aes = list("size" = 3), reverse = TRUE)
    ) +
    facet_grid(metag_complete ~ ., space = "free_y", scale = "free") +
    theme(
      strip.text.y = element_text(angle = 0)
    )
}

exp_design_plot(samp)
ggsave("exp_design_metagenomic.png", dpi = 700, width = 6.06, height = 6.62)
exp_design_plot(samp %>% filter(!no_interv))

## number of reads per sample
## reads <- data_frame(
##   "file" = fqs,
##   "n_reads" = NA
## )

## for (i in seq_along(fqs)) {
##   message(sprintf("Processing %s / %s", i, length(fqs)))
##   reads[i, "n_reads"]  <- system(sprintf("awk '{s++}END{print s/4}' %s", fqs[i]), intern = TRUE)
## }
## write.csv(reads, file = "../data/metagenomic/reads.csv", row.names = FALSE)
reads <- read.csv("../data/metagenomic/reads.csv") %>%
  mutate(
    Meas_ID = str_extract(file, "M[0-9]+"),
    forwards = grepl("1P", file)
  ) %>%
  left_join(meas) %>%
  left_join(samp) %>%
  filter(Samp_Type != "ExtrCont")

ggplot(reads) +
  geom_point(
    aes(
      x = Samp_Date, y = n_reads
    ),
    size = 0.5,
    alpha = 0.5
  ) +
  facet_wrap(~Subject, scale = "free_x")
ggsave("read_depths.png", dpi = 700)
