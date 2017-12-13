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
metag <- read_tsv("../data/metagenomic/merged/count_reads.tsv") %>%
  separate(species_id, c("Genus", "species", "strain_id"), remove = FALSE)

species_sums <- metag %>%
  select(starts_with("M")) %>%
  rowSums()
metag <- metag[species_sums != 0, ]


###############################################################################
## Functions used throughout
###############################################################################
count_taxa <- function(melted_abund, taxa_level = "Genus") {
  melted_abund %>%
    group_by_("sample", taxa_level) %>%
    summarise(total = sum(count)) %>%
    arrange(sample, desc(total)) %>%
    group_by(sample) %>%
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
  gather(sample, count, starts_with("M"))

melt_16s <- ps %>%
  get_taxa() %>%
  as.data.frame() %>%
  rownames_to_column("seq_id") %>%
  gather(sample, count, -seq_id) %>%
  left_join(taxa)

sums_metag <- count_taxa(melt_metag) %>%
  mutate(source = "metagenomic")
sums_16s <- count_taxa(melt_16s) %>%
  mutate(source = "16s")

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
## per sample plot of relative abundances between the two data sources
ggplot(genus) +
  geom_segment(
    aes(
      x = asinh(`16s`),
      y = sample,
      xend = asinh(metagenomic),
      yend = sample
    ),
    position = position_jitter(h = 0.5),
    alpha = 0.1,
    size = 0.2
  )

## barplot of relative abundances across genuses, from both data sets
keep_samples <- unique(sums_metag$sample)[1:10]
keep_genus <- levels(melt_genus$Genus)[1:100]
ggplot(
  melt_genus %>%
    filter(
      sample %in% keep_samples,
      Genus %in% keep_genus
    )
  ) +
  geom_bar(
    aes(x = Genus, y = rel_total, fill = source),
    stat = "identity",
    position = "dodge"
  ) +
  facet_grid(sample ~ .) +
  theme(axis.text.x = element_text(angle = 90, size = 6, hjust = 0))

## don't have any samples in common right now
intersect(sums_metag$sample %>% unique(), sums_16s$sample %>% unique())

## aggregate over all samples now
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

## scatterplot of logged relative abundances from the two sources
ggplot(
  ave_genus %>%
    spread(source, rel_total)
  ) +
  geom_abline(slope = 1, size = 0.1) +
  geom_text_repel(
    aes(x = `16s`, y = `metagenomic`, label = Genus),
    size = 2,
    force = 0.02
  ) +
  scale_x_log10() +
  scale_y_log10()
