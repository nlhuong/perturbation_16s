#!/usr/bin/env Rscript
##
## File description -------------------------------------------------------------
##
## Ordination plots for 16S data.
##
## author: nlhuong90@gmail.com
## date: 01/15/2018

#setwd("perturbation_16s/amplicon/")

###############################################################################
## Setup and data
###############################################################################

library(phyloseq)
library(dplyr)
library(tibble)
library(ggplot2)
library(viridis)
library(adaptiveGPCA)
library(doParallel)
library(foreach)
library(ggrepel)

options(stringsAsFactors = FALSE)
theme_set(theme_bw())
theme_update(text = element_text(size = 15))

path2data <- "../data/processed/"
path2figs <- "../figs"
path2out <- "../output"

ps0 <- readRDS(file.path(path2data, "perturb_physeq_filtered_27Dec.rds"))
remove_subjects <- c("AAA", "AAB", "AAN", "DAC")
ps <- subset_samples(ps0, !Subject %in% remove_subjects)

## Function to filter out taxa which occur in fewer than minPrev distinct sample 
## groups (e.g. Subjects).
filter_taxa <- function(physeq, minPrev = 1, group = "Subject") {
  smp <- data.frame(physeq@sam_data)
  seqtab <- as(physeq@otu_table, "matrix")
  group_prevalence <- sapply(1:ntaxa(physeq), function(i) {
    group_present <- smp[[group]][seqtab[i, ] > 0]
    length(unique(group_present))
  })
  physeq <- prune_taxa(taxa_names(physeq)[group_prevalence > minPrev],
                       physeq)
  return(physeq)
}


plot_projection <- function(data, xname, yname, labname = "Subject", size = 3,
                            color = NULL, eigs = NULL, ...){
  if(all(!is.null(color), color %in% colnames(data))){
    plt <- ggplot(data, aes_string(x = xname, y = yname, color = color)) 
  } else {
    plt <- ggplot(data, aes_string(x = xname, y = yname)) 
  }
  if(all(!is.null(labname), labname %in% colnames(data))) {
    if(size %in% colnames(data)) {
      plt <- plt + geom_text(aes_string(label = labname, size = size), ...)
    } else {
      plt <- plt + geom_text(aes_string(label = labname), size = size, ...)
    }
  } else {
    if(size %in% colnames(data)) {
      plt <- plt + geom_point(aes_string(size = size), ...)
    } else {
      plt <- plt + geom_point(size = size, ...)
    }
  }
  if(!is.null(eigs)){
    var.explained <- round(100 * eigs/sum(eigs), 2)
    plt <- plt + 
      labs(
        x = sprintf("Axis1 [%s%% variance]", var.explained[1]),
        y = sprintf("Axis2 [%s%% variance]", var.explained[2])
      ) +
      coord_fixed(sqrt(var.explained[2] / var.explained[1])) 
  }
  return(plt + theme(text = element_text(size = 20))) 
}

###############################################################################
## Standard PCA
###############################################################################

## Process data for each study arm separately 
## arms: ("Abx", "Diet_Abx", "CC_Diet", "CC_Diet_Abx", "NoIntv")

registerDoParallel(length(unique(ps@sam_data$Group)))
pca_study_arms <- foreach(i = seq_along(unique(ps@sam_data$Group))) %dopar% {
 g <- unique(ps@sam_data$Group)[i] 
 g_taxa <- taxa_names(ps)[ps@sam_data$Group == "NoIntv" | ps@sam_data$Group == g]
 ps_g <- prune_taxa(g_taxa, ps)
 ps_g <- filter_taxa(ps_g, minPrev = 2, group = "Subject")
 ps_g <- transform_sample_counts(ps_g, function(x) {asinh(x)})
 # PCA on centered data
 pca_g <- prcomp(scale(as(t(ps_g@otu_table), "matrix"), scale = FALSE))
 
 loadings <- pca_g$rotation[, 1:10] %>%
   as.data.frame() %>%
   rownames_to_column("Seq_ID") %>%
   left_join(data.frame(ps_g@tax_table) %>% rownames_to_column("Seq_ID"))
 
 scores <- pca_g$x[, 1:5] %>%
   as.data.frame() %>%
   rownames_to_column("Meas_ID") %>%
   left_join(ps_g@sam_data[, 1:30] %>% as.data.frame())%>%
   arrange(Group, Subject, Samp_Date)
 
 eig_df <- data.frame(
   eig_idx = 1:length(pca_g$sdev)
   ) %>%
   mutate(
     eig = pca_g$sdev^2,
     var_exp = 100*eig/sum(eig)
    )
 
 plt_scree <- ggplot(eig_df[1:10, ]) +
   geom_bar(aes(x = eig_idx, y = var_exp), stat = "identity")
 
 plt_scores <- plot_projection(
   scores, 
   xname = "PC1", yname = "PC2", 
   labname = "Subject", 
   size = 3, 
   color = "Subject", 
   eigs = eig_df$eig
 )
 
 plt_loadings <- plot_projection(
   loadings, 
   xname = "PC1", yname = "PC2", 
   size = "PC3",
   color = "Class", 
   eigs = eig_df$eig, 
   alpha = 0.3
 ) +
   geom_text_repel(
     data = loadings %>% 
       mutate(r2 = (PC1^2 + PC2^2)) %>%
       filter(r2 > quantile(r2, 0.95)), 
     aes(label = Genus)) + #, nudge_y = 0.001*max(sqrt(r2))) +
   theme(legend.position = "bottom")
 
  return(list(plt_scree, plt_scores, plt_loadings, pca_g))
}
names(pca_study_arms) <- unique(ps@sam_data$Group)

saveRDS(pca_study_arms, file = file.path(path2out, "pca_g.rds"))

pdf(file = file.path(path2figs, "study_arm_pca.pdf"))
for (g in names(pca_study_arms)) {
  g_pca <- pca_study_arms[[g]]
  plt <- gridExtra::grid.arrange(
    g_pca[[1]], g_pca[[2]], g_pca[[3]],
    ncol = 3, top = paste0("Study Arm: ", g))
  print(plt)
}
dev.off()

