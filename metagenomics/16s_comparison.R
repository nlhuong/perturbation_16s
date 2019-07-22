#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## The purpose of this script is to compare the abundance of species identified
## through metagenomic vs. 16s profiling. As is, this is a pretty coarse view,
## looking only at species assignments as estimated through MIDAS. A more
## careful analysis would consider the acutal counts of 16s sequences alone.
##
## author: sankaran.kris@gmail.com
## date: 12/11/2017

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

ps <- readRDS("../data/16s/perturb_physeq_8Dec.rds") %>%
  subset_samples(Subject %in% c("DBU", "DBV", "EAQ", "EAY", "EBF")) # samples with metagenomics
metag <- read_tsv("../data/metagenomic/merged/coverage.tsv") %>%
  separate(species_id, c("Genus", "species", "strain_id"), remove = FALSE)
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

species_sums <- metag %>%
  select(starts_with("M")) %>%
  rowSums()

bad_samples <- c("M3728", "M3673", "M3695", "M3204", "M3064", "M3109", "M3188",
                 "M3654")
metag <- metag[species_sums != 0, ]
metag <- metag[ !(colnames(metag) %in% bad_samples)]

###############################################################################
## Functions used throughout
###############################################################################
count_taxa <- function(melted_abund, taxa_level = "Genus") {
  melted_abund %>%
    group_by_("Meas_ID", taxa_level) %>%
    mutate(count = asinh(count)) %>%
    summarise(Samp_ID = Samp_ID[1], total = sum(count)) %>%
    arrange(Meas_ID, desc(total)) %>%
    group_by(Meas_ID) %>%
    mutate(rel_total = total / sum(total)) %>%
    ungroup()
}

###############################################################################
## Compare genus counts
###############################################################################
taxa <- tax_table(ps) %>%
  as.data.frame() %>%
  rownames_to_column("seq_id")

melt_metag <- metag %>%
  gather(Meas_ID, count, starts_with("M")) %>%
  left_join(meas %>% select("Meas_ID", "Samp_ID"))

melt_16s <- ps %>%
  get_taxa() %>%
  as.data.frame() %>%
  rownames_to_column("seq_id") %>%
  gather(Meas_ID, count, -seq_id) %>%
  left_join(taxa %>% select("seq_id", "Genus")) %>%
  left_join(meas %>% select("Meas_ID", "Samp_ID"))

sums_metag <- count_taxa(melt_metag) %>%
  mutate(source = "metagenomic")
sums_16s <- count_taxa(melt_16s) %>%
  mutate(source = "16s")
rm(melt_16s)
rm(melt_metag)

genus_levels <- c(sums_metag$Genus, levels(sums_16s$Genus)) %>%
  unique()
melt_genus <- bind_rows(sums_metag, sums_16s) %>%
  mutate(Genus = factor(Genus, genus_levels))

genus <- melt_genus %>%
  spread(source, rel_total) %>%
  mutate(
    `16s` = ifelse(is.na(`16s`), 0, `16s`),
  ) %>%
  filter(!is.na(metagenomic))

###############################################################################
## Plot genus count comparions
###############################################################################
## barplot of relative abundances across genuses, from both data sets
keep_ids <- meas %>%
  filter(Meas_ID %in% c(colnames(metag), sample_names(ps))) %>%
  group_by(Samp_ID) %>%
  mutate(n_types = n()) %>%
  filter(n_types > 1) %>%
  arrange(Samp_ID) %>%
  .[["Meas_ID"]]

keep_genus <- levels(melt_genus$Genus)[1:100]
ggplot(
  melt_genus %>%
    filter(
      Meas_ID %in% keep_ids[1:10],
      Genus %in% keep_genus
    )
  ) +
  geom_bar(
    aes(x = Genus, y = rel_total, fill = source),
    stat = "identity",
    position = "dodge"
  ) +
  facet_grid(Samp_ID ~ ., scale = "free_y") +
  theme(axis.text.x = element_text(angle = 90, size = 6, hjust = 0))

## scatterplot of abundances from the two sources, per sample
scatter_data <- melt_genus %>%
  filter(Genus %in% keep_genus) %>%
  left_join(samp) %>%
  select(-Meas_ID, -total) %>%
  spread(source, rel_total)

genus_centroid <- scatter_data %>%
  group_by(Subject, Genus, Abx_Interval) %>%
  summarise(
    `16s` = median(`16s`),
    metagenomic = median(metagenomic, na.rm = TRUE)
  )

ggplot() +
  geom_abline(alpha = 0.4, size = 0.2) +
  geom_point(
    data = scatter_data,
    aes(x = log(`16s`), y = log(`metagenomic`), col = Abx_Interval),
    size = 0.5, alpha = 0.2
  ) +
  geom_point(
    data = genus_centroid,
    aes(x = log(`16s`), y = log(`metagenomic`), fill = Abx_Interval),
    size = 1, pch = 21, stroke = 0.2
  ) +
  geom_text_repel(
    data = genus_centroid %>%
      filter(
        `16s` > 0,
        metagenomic > 0,
        abs(log(`16s`) - log(metagenomic)) > 3
      ),
    aes(x = log(`16s`), y = log(`metagenomic`), label = Genus),
    size = 2,
    force = 0.2
  ) +
  facet_wrap(~Subject) +
  theme(legend.position = "bottom")
ggsave("comparison_scatter.png", width = 5.28, height = 3.54, dpi = 600)
