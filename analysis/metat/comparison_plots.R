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

###############################################################################
## Data Loading
###############################################################################

# loading sample info
meas <- read_xlsx("data/Mapping_Files_22Jan2018.xlsx", "Meas", skip = 1) %>%
  rename(Samp_ID = SampID)
interv_levs <- c("NoInterv", "PreDiet", "MidDiet", "PostDiet", "PreCC", "MidCC",
                 "PostCC", "PreAbx", "MidAbx", "PostAbx")
samp <- read_xlsx("data/Mapping_Files_22Jan2018.xlsx", "Samp", skip = 1) %>%
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
# MetatT Data
metat <- read_csv("data/metatranscriptomic/counts/Annotated_Relman_RNAseq_16_refseq_genes.csv") 
geneSums <- metat %>% select(starts_with("M")) %>% rowSums()
metat <- metat[geneSums > 10, ]
RPKM <- sweep(metat%>%select(starts_with("M")), 2, colSums(metat%>%select(starts_with("M"))), "/")
RPKM <- 10^9 * RPKM
RPKM <- sweep(RPKM, 1, metat$Length, "/")
rownames(RPKM) <- metat$GeneID
RPKM <- RPKM[rowSums(RPKM) > 5, ]
metatInfo <- metat %>% filter(GeneID %in% rownames(RPKM)) %>% select(GeneID, Length, Organism, Function)
rm(metat)

metatInfo <- metatInfo %>%
  separate(Organism, c("Genus", "species"), sep = " ", remove = FALSE,
           extra = "merge", fill = "right") %>%
  mutate(species_id = Organism)

smp_names <- colnames(RPKM) 
colnames(RPKM) <- gsub("\\_(.*)", "", smp_names)
metat <- cbind(metatInfo, RPKM)
rm(RPKM, metatInfo)

# 16S Data
ps <- readRDS("data/16s/perturb_physeq_22Jan18.rds") 

# MetaG Data
metag <- read_tsv("data/metagenomic/coverage.tsv") %>%
  separate(species_id, c("Genus", "species", "strain_id"), 
           remove = FALSE, extra = "merge")
bad_samples <- c("M3728", "M3673", "M3695", "M3204", "M3064", "M3109", "M3188",
                 "M3654")
metag <- metag[ !(colnames(metag) %in% bad_samples)]


# samples with in all 3 datasets:
keep_samples <- meas %>%
  filter(Meas_ID %in% c(colnames(metag), colnames(metat_RPKM), sample_names(ps))) %>%
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
metag <- metag[metag_rowsums > 0, ]

metat_rowsums <- metat %>% 
  select(starts_with("M")) %>%
  rowSums()
metat <- metat[metat_rowsums > 1, ]

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

sums_metat <- count_taxa(melt_metat) %>%
  mutate(source = "metatranscriptomic")
rm(melt_metat)

melt_metag <- metag %>%
  gather(Meas_ID, count, starts_with("M")) %>%
  left_join(meas %>% select("Meas_ID", "Samp_ID"))
sums_metag <- count_taxa(melt_metag) %>%
  mutate(source = "metagenomic")
rm(melt_metag)

melt_16s <- ps %>%
  get_taxa() %>%
  as.data.frame() %>%
  rownames_to_column("seq_id") %>%
  gather(Meas_ID, count, -seq_id) %>%
  left_join(taxa %>% select("seq_id", "Genus")) %>%
  left_join(meas %>% select("Meas_ID", "Samp_ID"))

sums_16s <- count_taxa(melt_16s) %>%
  mutate(source = "16s")
rm(melt_16s)


# melt_genus <- bind_rows(sums_metat, sums_metag, sums_16s) 
# 

melt_genus <- bind_rows(sums_metat, sums_16s, sums_metag)
genus <- melt_genus  %>%
  select(-Meas_ID, -total) %>%
  spread(key = source, value = rel_total) %>%
  mutate(
    `16s` = ifelse(is.na(`16s`), 0, `16s`)
  ) 

genus_levels <- genus %>%
  select(-Samp_ID) %>%
  group_by(Genus) %>%
  summarise_all(mean) %>%
  mutate(mean_dataset = `16s` + metatranscriptomic + metagenomic) %>%
  arrange(desc(mean_dataset))

## scatterplot of abundances from the two sources, per sample
scatter_data <- melt_genus  %>%
  filter(Genus %in% genus_levels$Genus) %>%
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

melt_genus <- melt_genus %>%
  mutate(Genus = factor(Genus, genus_levels$Genus))



###############################################################################
## Plot genus count comparions
###############################################################################


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


plot_comparison <- function(scatter_data, genus_centroid,
                            x = "16s", y = "metatranscriptomic", 
                            color = "Diet_Interval", label = "Genus") {
  
  genus_centroid <- genus_centroid %>% as.data.frame() 
  scatter_data <- scatter_data %>% as.data.frame()
  
  genus_centroid$x <- genus_centroid[[x]]
  genus_centroid$y <- genus_centroid[[y]]
  scatter_data$x <- scatter_data[[x]]
  scatter_data$y <- scatter_data[[y]]
  
  genus_centroid <- genus_centroid %>%
    filter(!is.na(x) & !is.na(y))
  scatter_data <- scatter_data %>%
    filter(!is.na(x) & !is.na(y))
  
  data_repel <-  genus_centroid %>%
    filter(
      x > 0, y > 0, 
      abs(log(x)- log(y)) > 2.5
    ) %>%
    mutate(x = log(x), y = log(y) )
  
  genus_centroid[, c("x", "y")] <- log(genus_centroid[, c("x", "y")])
  scatter_data[, c("x", "y")] <- log(scatter_data[, c("x", "y")])
  
  ggplot() +
    geom_abline(alpha = 0.4, size = 0.2) +
    geom_point(
      data = scatter_data,
      aes_string("x", "y", col = color),
      size = 0.9, alpha = 0.5
    ) +
    geom_point(
      data = genus_centroid,
      aes_string(x = "x", y = "y", fill = color),
      size = 1.5, pch = 21, stroke = 0.2
    ) +
    geom_text_repel(
      data = data_repel,
      aes_string(x = "x", y = "y", label = label),
      size = 3,
      force = 0.2
    ) +
    facet_wrap(~Subject) + theme_scatter + #xlim(log(-0.5), NA) +
    xlab(paste0("log(", x, ")")) + ylab(paste0("log(", y, ")"))
}


pdf("figs/comparison_scatter.pdf", width = 8, height = 6)
plot_comparison(scatter_data, genus_centroid,
                x = "16s", y = "metatranscriptomic", 
                color = "Diet_Interval", label = "Genus") 
plot_comparison(scatter_data, genus_centroid,
                x = "metagenomic", y = "metatranscriptomic", 
                color = "Diet_Interval", label = "Genus") 
plot_comparison(scatter_data, genus_centroid,
                x = "16s", y = "metagenomic", 
                color = "Diet_Interval", label = "Genus") 
dev.off()
