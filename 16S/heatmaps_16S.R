#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Heatmaps for 16S data.
##
## author: nlhuong90@gmail.com
## date: 01/03/2018

###############################################################################
## Setup and data
###############################################################################
#setwd("perturbation_16s/amplicon/")

## Load libraries
.packages <- c("gtable","grid", "ggplot2", "viridis", "RColorBrewer", 
               "phyloseq", "tidyverse", "pheatmap", "reshape2")
sapply(.packages, require, character.only = TRUE)

options(stringsAsFactors = FALSE)

## Set custom plot defaults
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
  legend.text = element_text(size = 9),
  axis.text = element_text(size = 9),
  axis.title = element_text(size = 8),
  strip.background = element_blank(),
  strip.text = element_text(size = 9),
  legend.key = element_blank()
)

## Read data
path_to_data <- "../data/processed/perturb_physeq_filtered_27Dec.rds"
ps <- readRDS(path_to_data)

## Plotting functions 
plot_heatmap <- function(mabund, fill_type = NULL, taxrank = NULL, 
                         ntax = 20, na.tax = c("keep", "remove", "select")) {
  if(!is.null(taxrank)) {
    mabund$taxrank <- mabund[, taxrank]
    tax_levels <- mabund %>%
      group_by(taxrank) %>%
      summarise(prev = sum(abundance > 0, na.rm =TRUE)) %>%
      arrange(desc(prev)) %>%
      .[["taxrank"]]
    tax_levels <- tax_levels[!is.na(tax_levels)]
    mabund <- mutate(mabund, taxrank = factor(taxrank, tax_levels[1:ntax]))
    if(na.tax == "remove") 
      mabund <- filter(mabund, !is.na(taxrank))
    if(na.tax == "select") 
      mabund <- filter(mabund, is.na(taxrank))
  }
  
  if(is.null(fill_type)){
    plt <- ggplot(mabund) +
      geom_raster( 
        aes_string(x = "Seq_ID", y = "Meas_ID", fill = "abundance")
      ) +
      scale_fill_viridis()
  } else {
    plt <- ggplot(mabund)  + 
      geom_raster( 
        aes_string(x = "Seq_ID", y = "Meas_ID", 
                   alpha = "abundance", fill = fill_type)
      ) +
      scale_alpha(range = c(0, 1)) +
      scale_fill_discrete()
  }
  if(!is.null(taxrank)) {
    plt <- plt + 
      facet_grid(Group + Subject ~ taxrank, switch = "y",
                 scales = "free", space = "free")
  } else {
    plt <- plt + 
      facet_grid(Group +Subject ~ ., switch = "y",
                 scales = "free_y", space = "free")
  }
  
  plt + theme(
    axis.text = element_blank(),
    axis.title.y = element_blank(),
    strip.text.x = element_text(angle = 90, hjust = 0),
    strip.text.y = element_text(angle = 90),
    legend.position = "bottom"
  )
    
}

###############################################################################
## Prepare data for heatmaps 
###############################################################################

## Filter samples
remove_subjects <- c("AAA", "AAB", "AAN", "DAC") # not in longitudinal study
ps <- subset_samples(ps, !Subject %in% remove_subjects)

## Filter taxa
minTaxaSubjPrev <- 10
num_subj_present <- function(physeq, subject_indicator = "Subject") {
  smp <- physeq@sam_data %>%
    as.data.frame()
  seqtab <- as(physeq@otu_table, "matrix")
  sapply(1:ntaxa(physeq), function(i) {
    subj.present <- smp[[subject_indicator]][seqtab[i, ] > 0]
    length(unique(subj.present))})
}
subjects_prevalence <- num_subj_present(ps)
ps1 <- subset_taxa(ps, subjects_prevalence >= minTaxaSubjPrev) # 1036 x 2932

## Process genus and species names for plotting
sub_na <- function(x, taxtab, subs){
  taxtab <- data.frame(taxtab)
  res <- c()
  for (sub in subs) {
    res <- sapply(1:nrow(taxtab), function(i)
      ifelse(!is.na(x[i]), x[i], taxtab[i, sub]))
  }
  return(res)
}

## Trandform (asinh) data
psIHS <- transform_sample_counts(ps1, function(x) {log(x + sqrt(x^2+1))}) 

## Combine assignment from RDP and Silva reference db
taxtab <- psIHS@tax_table %>%
  as.data.frame() %>%
  rownames_to_column("Seq_ID")
taxtab <- taxtab %>%
  mutate(
    genus = sub_na(taxtab$Genus, taxtab, "GenusSilva"),
    species = sub_na(taxtab$Species, taxtab, c("SpeciesSilva", "Seq_ID"))
  )

## Convert to long data format
mabund <- psIHS@otu_table %>%
  melt(value.name = "abundance", varnames  = c("Seq_ID", "Meas_ID")) %>%
  left_join(psIHS@sam_data %>% select(Meas_ID:Samp_Date, SampID)) %>%
  left_join(taxtab %>% select(Seq_ID, genus, species)) 

## Set sequence ordering by pheatmap (hclust)
abundheat <- pheatmap(
  as(psIHS@otu_table, "matrix"),
  cluster_cols = FALSE,
  silent = TRUE)
seq_order <- abundheat$tree_row$order
mabund <- mabund %>%
  mutate(Seq_ID = factor(Seq_ID, rownames(psIHS@otu_table)[seq_order]))

## Set measurement ordering levels
meas_levels <- psIHS@sam_data %>%
  as.data.frame() %>%
  select(Meas_ID, Subject, Samp_Date) %>%
  arrange(Subject, desc(Samp_Date)) %>%
  .[["Meas_ID"]]
mabund <- mutate(mabund, Meas_ID = factor(Meas_ID, meas_levels))


###############################################################################
## Generate heatmaps
###############################################################################

pdf(file = "../figs/heatmaps_16S_asinh.pdf", width=12, height=14)
for (g in unique(mabund$Group)) {
  g.mabund <- mabund %>% filter(Group == g)
  pheat <- plot_heatmap(g.mabund,
    fill_type = NULL, taxrank = NULL) +
    ggtitle(paste("Group:", g))
  print(pheat)
}
dev.off()




