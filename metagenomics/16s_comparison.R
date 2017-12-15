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

ps <- readRDS("../data/16s/perturb_physeq_8Dec.rds")
metag <- read_tsv("../data/metagenomic/merged/coverage.tsv") %>%
  separate(species_id, c("Genus", "species", "strain_id"), remove = FALSE)
meas <- read_xlsx("../data/Mapping_Files_7bDec2017.xlsx", "Meas", skip = 1)

species_sums <- metag %>%
  select(starts_with("M")) %>%
  rowSums()
metag <- metag[species_sums != 0, ]

###############################################################################
## Functions used throughout
###############################################################################
count_taxa <- function(melted_abund, taxa_level = "Genus") {
  melted_abund %>%
    group_by_("Meas_ID", taxa_level) %>%
    summarise(SampID = SampID[1], total = sum(count)) %>%
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
  left_join(meas %>% select("Meas_ID", "SampID"))

melt_16s <- ps %>%
  get_taxa() %>%
  as.data.frame() %>%
  rownames_to_column("seq_id") %>%
  gather(Meas_ID, count, -seq_id) %>%
  left_join(taxa %>% select("seq_id", "Genus")) %>%
  left_join(meas %>% select("Meas_ID", "SampID"))

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
  group_by(SampID) %>%
  mutate(n_types = n()) %>%
  filter(n_types > 1) %>%
  arrange(SampID) %>%
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
  facet_grid(SampID ~ ., scale = "free_y") +
  theme(axis.text.x = element_text(angle = 90, size = 6, hjust = 0))

## scatterplot of abundances from the two sources, per sample
ggplot(
  melt_genus %>%
    filter(
      Meas_ID %in% keep_ids[1:12],
      Genus %in% keep_genus
    ) %>%
    select(-Meas_ID, -total) %>%
    spread(source, rel_total) %>%
    filter(`16s` > 1e-3 | `metagenomic` > 1e-3)
  ) +
  geom_vline(xintercept = 0, alpha = 0.4, size = 0.2) +
  geom_hline(yintercept = 0, alpha = 0.4, size = 0.2) +
  geom_abline(alpha = 0.4, size = 0.2) +
  geom_text_repel(
    aes(x = sqrt(`16s`), y = sqrt(`metagenomic`), label = Genus),
    size = 2,
    force = 0.02
  ) +
  facet_wrap(~SampID, scale = "free")

## aggregate over all SampIDs now
ave_genus <- melt_genus %>%
    filter(Genus %in% keep_genus) %>%
    group_by(Genus, source) %>%
    summarise(rel_total = mean(rel_total, na.rm = TRUE))

ggplot(ave_genus) +
  geom_bar(
    aes(x = Genus, y = rel_total, fill = source),
    stat = "identity",
    position = "dodge"
  ) +
  theme(axis.text.x = element_text(angle = 90, size = 6, hjust = 0))

ggplot(
  ave_genus %>%
    filter(Genus %in% keep_genus) %>%
    spread(source, rel_total) %>%
    filter(`16s` > 1e-3 | `metagenomic` > 1e-3)
  ) +
  geom_vline(xintercept = 0, alpha = 0.4, size = 0.2) +
  geom_hline(yintercept = 0, alpha = 0.4, size = 0.2) +
  geom_abline(alpha = 0.4, size = 0.2) +
  geom_text_repel(
    aes(x = sqrt(`16s`), y = sqrt(`metagenomic`), label = Genus),
    size = 2,
    force = 0.02
  )
