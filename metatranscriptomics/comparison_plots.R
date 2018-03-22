#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## The purpose of this script is to compare the abundance of species identified
## through metatranscriptomics vs. metagenomic vs. 16s profiling. 
## We compare the counts generated with metratranscriptomics pipeline with
## the ones from metagenomics data processing with MIDAS and also with 
## the amplicon 16s counts processed with DADA2.
##
## author: nlhuong90@gmail.com
## date: 03/15/2018

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
theme_updates <- theme(
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

# loading sample info
meas <- read_xlsx("../data/Mapping_Files_22Jan2018.xlsx", "Meas", skip = 1) %>%
  rename(Samp_ID = SampID)
interv_levs <- c("NoInterv", "PreDiet", "MidDiet", "PostDiet", "PreCC", "MidCC",
                 "PostCC", "PreAbx", "MidAbx", "PostAbx")
samp <- read_xlsx("../data/Mapping_Files_22Jan2018.xlsx", "Samp", skip = 1) %>%
  mutate(
    Diet_Interval = ifelse(Diet_Interval == "NA", "NoInterv", Diet_Interval),
    CC_Interval = ifelse(CC_Interval == "NA", "NoInterv", CC_Interval),
    Abx_Interval = ifelse(Abx_Interval == "NA", "NoInterv", Abx_Interval),
    Diet_Interval = factor(Diet_Interval, interv_levs),
    CC_Interval = factor(CC_Interval, interv_levs),
    Abx_Interval = factor(Abx_Interval, interv_levs)
  ) %>%
  filter(Samp_Type != "ExtrCont")

# loading counts
ps <- readRDS("../data/16s/perturb_physeq_22Jan18.rds") 
metag <- read_tsv("../data/metagenomic/coverage.tsv") %>%
  separate(species_id, c("Genus", "species", "strain_id"), 
           remove = FALSE, extra = "merge")
bad_samples <- c("M3728", "M3673", "M3695", "M3204", "M3064", "M3109", "M3188",
                 "M3654")
metag <- metag[ !(colnames(metag) %in% bad_samples)]
metat <- read_csv("../data/metatranscriptomic/DBUr_Sub_refseq_organism.csv") %>%
  select(-X1) %>%
  separate(feature_name, c("Genus", "species"), sep = " ", remove = FALSE,
           extra = "merge", fill = "right")
colnames(metat) <- gsub("\\_(.*)", "", colnames(metat))
metat <- metat  %>%  rename(species_id = feature)

# samples with in all 3 datasets:
keep_samples <- meas %>%
  filter(Meas_ID %in% c(colnames(metag), colnames(metat), sample_names(ps))) %>%
  group_by(Samp_ID) %>%
  mutate(n_types = n()) %>%
  filter(n_types == 3) %>%
  arrange(Samp_ID) 

length(unique(keep_samples$Samp_ID))

keep_meas <- keep_samples %>%
  .[["Meas_ID"]]

metat <- metat[ , c("species_id", "Genus", "species", intersect(colnames(metat), keep_meas))]
metag <- metag[, c("species_id", "Genus", "species", intersect(colnames(metag), keep_meas))]
ps <- subset_samples(ps, sample_names(ps) %in% keep_meas)


# filter taxa
metag_rowsums <- metag %>%
  select(starts_with("M")) %>%
  rowSums()
metag <- metag[metag_rowsums != 0, ]

metat_rowsums <- metat %>% 
  select(starts_with("M")) %>%
  rowSums()
metat <- metat[metat_rowsums != 0, ]

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

melt_metat <- metat %>%
  gather(Meas_ID, count, starts_with("M")) %>%
  left_join(meas %>% select("Meas_ID", "Samp_ID"))

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

sums_metat <- count_taxa(melt_metat) %>%
  mutate(source = "metatranscriptomic")

sums_metag <- count_taxa(melt_metag) %>%
  mutate(source = "metagenomic")

sums_16s <- count_taxa(melt_16s) %>%
  mutate(source = "16s")

rm(melt_16s, melt_metag, melt_metat)

melt_genus <- bind_rows(sums_metat, sums_metag, sums_16s) 

genus <- melt_genus %>%
  select(-Meas_ID, -total) %>%
  spread(key = source, value = rel_total) %>%
  mutate(
    `16s` = ifelse(is.na(`16s`), 0, `16s`)
  ) %>%
  filter(!is.na(metagenomic), !is.na(metatranscriptomic))

genus_levels <- genus %>%
  select(-Samp_ID) %>%
  group_by(Genus) %>%
  summarise_all(mean) %>%
  mutate(mean_dataset = `16s` + metagenomic + metatranscriptomic) %>%
  arrange(desc(mean_dataset))

melt_genus <- melt_genus %>%
  mutate(Genus = factor(Genus, genus_levels$Genus)) 


###############################################################################
## Plot genus count comparions
###############################################################################
keep_genus <- genus_levels$Genus
## scatterplot of abundances from the two sources, per sample
scatter_data <- melt_genus %>%
  filter(Genus %in% keep_genus) %>%
  left_join(samp) %>%
  select(-Meas_ID, -total) %>%
  spread(source, rel_total)

genus_centroid <- scatter_data %>%
  group_by(Subject, Genus, Diet_Interval) %>%
  summarise(
    `16s` = median(`16s`),
    metagenomic = median(metagenomic, na.rm = TRUE),
    metatranscriptomic = median(metatranscriptomic, na.rm = TRUE)
  )


pdf("../figs/comparison_scatter.pdf", width = 8, height = 6)

ggplot(
  melt_genus %>%
    filter(
      Samp_ID %in% keep_samples$Samp_ID[1:10],
      !is.na(Genus),
      Genus %in% keep_genus
    )
) +
  geom_bar(
    aes(x = Genus, y = rel_total, fill = source),
    stat = "identity",
    position = "dodge"
  ) +
  facet_grid(Samp_ID ~ ., scale = "free_y") +
  theme(axis.text.x = element_text(angle = 90, size = 6, hjust = 0)) +
  theme_updates 
  

theme_scatter <- theme(
  legend.position = "bottom",  
  text = element_text(size = 12),
  strip.text = element_text(size = 11),
  strip.background = element_blank()
)

ggplot() +
  geom_abline(alpha = 0.4, size = 0.2) +
  geom_point(
    data = scatter_data,
    aes(x = log(`16s`), y = log(`metagenomic`), col = Diet_Interval),
    size = 0.9, alpha = 0.5
  ) +
  geom_point(
    data = genus_centroid,
    aes(x = log(`16s`), y = log(`metagenomic`), fill = Diet_Interval),
    size = 1.5, pch = 21, stroke = 0.2
  ) +
  geom_text_repel(
    data = genus_centroid %>%
      filter(
        `16s` > 0,
        metagenomic > 0,
        abs(log(`16s`) - log(metagenomic)) > 2.5
      ),
    aes(x = log(`16s`), y = log(`metagenomic`), label = Genus),
    size = 3,
    force = 0.2
  ) +
  facet_wrap(~Subject) + theme_scatter


ggplot() +
  geom_abline(alpha = 0.4, size = 0.2) +
  geom_point(
    data = scatter_data,
    aes(x = log(`16s`), y = log(`metatranscriptomic`), col = Diet_Interval),
    size = 0.9, alpha = 0.5
  ) +
  geom_point(
    data = genus_centroid,
    aes(x = log(`16s`), y = log(`metatranscriptomic`), fill = Diet_Interval),
    size = 1.5, pch = 21, stroke = 0.2
  ) +
  geom_text_repel(
    data = genus_centroid %>%
      filter(
        `16s` > 0,
        metagenomic > 0,
        abs(log(`16s`) - log(metatranscriptomic)) > 2.5
      ),
    aes(x = log(`16s`), y = log(`metatranscriptomic`), label = Genus),
    size = 3,
    force = 0.2
  ) +
  facet_wrap(~Subject) +
theme_scatter

ggplot() +
  geom_abline(alpha = 0.4, size = 0.2) +
  geom_point(
    data = scatter_data,
    aes(x = log(`metagenomic`), y = log(`metatranscriptomic`), col = Diet_Interval),
    size = 0.9, alpha = 0.5
  ) +
  geom_point(
    data = genus_centroid,
    aes(x = log(`metagenomic`), y = log(`metatranscriptomic`), fill = Diet_Interval),
    size = 1.5, pch = 21, stroke = 0.2
  ) +
  geom_text_repel(
    data = genus_centroid %>%
      filter(
        `16s` > 0,
        metagenomic > 0,
        abs(log(`16s`) - log(metatranscriptomic)) > 2.5
      ),
    aes(x = log(`16s`), y = log(`metatranscriptomic`), label = Genus),
    size = 3,
    force = 0.2
  ) +
  facet_wrap(~Subject) + theme_scatter

dev.off()