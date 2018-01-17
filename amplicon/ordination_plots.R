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
rm(list = ls())

library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(viridis)
library(phyloseq)
library(adaptiveGPCA)
library(doParallel)
library(foreach)

rm(list = ls())
options(stringsAsFactors = FALSE)
theme_set(theme_bw())
theme_update(text = element_text(size = 15))

path2data <- "../data/processed/"
path2figs <- "../figs"
path2out <- "../output"

## Load data
ps0 <- readRDS(file.path(path2data, "perturb_physeq_filtered_27Dec.rds"))
remove_subjects <- c("AAA", "AAB", "AAN", "DAC")
ps <- subset_samples(ps0, !Subject %in% remove_subjects)
ps@sam_data <- ps@sam_data[, 1:30]


## Function to filter out taxa which occur in fewer than minPrev distinct sample 
## groups (e.g. Subjects).
filter_taxa <- function(physeq, minPrev = 1, group = "Subject") {
  if(length(unique(physeq@sam_data[[group]])) == 1) return(physeq)
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


get_ordinations <- function(physeq, group = NULL, method = "pca", ncores = 2) {
  if(is.null(group) || !group %in% colnames(physeq@sam_data))
    physeq@sam_data[[group]] <- "group1"
  groups <- unique(physeq@sam_data[[group]])
  
  registerDoParallel(min(ncores, length(groups))) 
  res <- foreach(i = seq_along(groups)) %dopar% {
    g <- unique(physeq@sam_data[[group]])[i] 
    g_samples <- rownames(physeq@sam_data)[physeq@sam_data[[group]] == g] 
    g_physeq <- prune_samples(g_samples, physeq)
    g_physeq <- filter_taxa(g_physeq, minPrev = 1, group = "Subject")

    g_physeq <- transform_sample_counts(g_physeq, function(x) {asinh(x)})
    g_taxtab <- data.frame(g_physeq@tax_table) %>% rownames_to_column("Seq_ID")
    
    eigs <- loadings <- scores <- NULL
    if(tolower(method) == "pca") {
      g_ord <- prcomp(scale(as(t(g_physeq@otu_table), "matrix"), scale = FALSE))
      scores <- as.data.frame(g_ord$x)
      loadings <- as.data.frame(g_ord$rotation)
      eigs <- data.frame(
        eig_idx = 1:length(g_ord$sdev), 
        eig = g_ord$sdev^2) %>%
        mutate(var_exp = 100*eig/sum(eig))
    } else if (tolower(method) == "agpca") {
      pp <- processPhyloseq(g_physeq)
      g_ord <- adaptivegpca(pp$X, pp$Q, k = 3)
      scores <- as.data.frame(g_ord$U)
      loadings <- as.data.frame(g_ord$QV)
      eigs <- data.frame(
        eig_idx = 1:length(g_ord$vars),
        eig = g_ord$vars, vars_exp = g_ord$vars)
    } else if (tolower(method) == "tsne") {
      brayD <- phyloseq::distance(g_physeq, method = "bray")
      g_ord <- Rtsne::Rtsne(brayD, is_distance = TRUE, dims = 2, pca = TRUE,
                   eta = 1, exaggeration_factor = nsamples(g_physeq)/10)
      scores <- as.data.frame(g_ord$Y)
    } else {
      stop("method not supported")
    }
    if(!is.null(loadings)){
      colnames(loadings) <- paste0("Axis", 1:ncol(loadings))
      loadings <- loadings %>%
        rownames_to_column("Seq_ID") %>%
        left_join(g_taxtab)
    }
    colnames(scores) <- paste0("Axis", 1:ncol(scores))
    scores <- scores %>%
      rownames_to_column("Meas_ID") %>%
      left_join(g_physeq@sam_data %>% as.data.frame())%>%
      arrange(Group, Subject, Samp_Date, DayFromStart, Timeline)

    return(list(eigs = eigs, scores = scores, loadings = loadings))
  }
  names(res) <- groups
  return(res)
}


plot_projection <- function(data, xname, yname, labname = NULL, size = 3,
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

plot_scores <- function(scores, size = 3, eigs = NULL) {
  plot_projection(
    scores, 
    xname = "Axis1", 
    yname = "Axis2", 
    labname = "Subject",
    size = size,
    color = "Subject", 
    eigs = eigs)
} 

plot_scores_time <- function(scores, size = 3, eigs = NULL){
  suppressWarnings(
    plot_projection(
      scores, 
      xname = "Axis1", 
      yname = "Axis2", 
      color = "DayFromStart", 
      size = size,
      eigs = eigs) + 
      geom_point(
        data = scores %>% filter(Timeline != "typical"), 
        aes(fill = Timeline, label = Subject),
        color = "black",
        size = 1.5*size,
        shape = 23,
        lwd=2) +
      scale_color_viridis() +
      scale_fill_brewer(palette = "Oranges")
  )
} 

plot_loadings <- function(loadings, size = 3, eigs = NULL, 
                          color = "Class", label = "Genus"){
  
  loadings_fltr <- loadings %>%
    mutate(r2 = (Axis1^2 + Axis2^2)) %>%
    filter(
      r2 > quantile(r2, 0.97),
      Axis1 > quantile(Axis1, 0.99) | Axis1 < quantile(Axis1, 0.01),
      Axis2 > quantile(Axis2, 0.99) | Axis2 < quantile(Axis2, 0.01))
    plot_projection(
      loadings, 
      xname = "Axis1", 
      yname = "Axis2", 
      size = size,
      color = "Class", 
      eigs = eigs, 
      alpha = 0.3
    ) +
    geom_text_repel(
        data = loadings_fltr,
        aes_string(label = label) 
  )
}

generate_pdf <- function(ord_lst, filename, width = 15, height = 30) {
  ord_plts <- lapply(ord_lst, function(g_ord) {
    eigs <- p_scree <- p_scores<- NULL
    if(!is.null(g_ord$eigs)){
      eigs <- g_ord$eigs$eig
      p_scree <- ggplot(g_ord$eigs[1:10, ]) +
        geom_bar(aes(x = eig_idx, y = var_exp), stat = "identity")
    }
    if(length(unique(g_ord$scores[["Subject"]])) > 1) 
      p_scores <- plot_scores(g_ord$scores, eigs = eigs)
    p_scores_time <- plot_scores_time(g_ord$scores, eigs = eigs)
    p_loadings <- plot_loadings(g_ord$loadings, eigs = eigs)
    return(list(scree = p_scree, loadings = p_loadings, 
                scores = p_scores, scores_time = p_scores_time))
  })
  names(ord_plts) <- names(ord_lst)
  
  pdf(file = filename, width = width, height = height)
  for (g in names(ord_plts)) {
    plt <- ord_plts[[g]]
    plt <- plt[!sapply(plt, is.null)]
    grid.arrange(grobs = plt, ncol = 1, top = paste0("Study Arm: ", g))
  }
  dev.off()
}

###############################################################################
## Standard PCA
###############################################################################

## Process data for each study arm separately 
## arms: ("Abx", "Diet_Abx", "CC_Diet", "CC_Diet_Abx", "NoIntv")
pca_res <- get_ordinations(ps, group = "Group", method = "pca")
names(pca_res) <- unique(ps@sam_data[["Group"]])
saveRDS(pca_res, file = file.path(path2out, "study_arm_pca.rds"))

generate_pdf(ord_lst = pca_res, 
             filename = file.path(path2figs, "study_arm_pca.pdf"), 
             width = 15, height = 30)

## Process each subject separately
pca_subject_res <- get_ordinations(ps, group = "Subject", method = "PCA")
saveRDS(pca_subject_res, file = file.path(path2out, "subject_pca.rds"))

generate_pdf(ord_lst = pca_subject_res, 
             filename = file.path(path2figs, "subject_pca.pdf"), 
             width = 15, height = 25)

###############################################################################
## Adaptive GPCA
###############################################################################

## Process data for each study arm separately 
agpca_res <- get_ordinations(ps, group = "Group", method = "agpca")
saveRDS(agpca_res, file = file.path(path2out, "study_arm_agpca.rds"))

generate_pdf(ord_lst = agpca_res, 
             filename = file.path(path2figs, "study_arm_agpca.pdf"), 
             width = 15, height = 30)

## Process each subject separately
agpca_subject_res <- get_ordinations(ps, group = "Subject", method = "agpca")
saveRDS(agpca_subject_res, file = file.path(path2out, "subject_pca.rds"))

generate_pdf(ord_lst = pca_subject_res, 
             filename = file.path(path2figs, "subject_pca.pdf"), 
             width = 15, height = 25)

###############################################################################
## tSNE
###############################################################################

## Process data for each study arm separately 
tsne_res <- get_ordinations(ps, group = "Group", method = "tsne")
saveRDS(tsne_res, file = file.path(path2out, "study_arm_tsne.rds"))

generate_pdf(ord_lst = tsne_res, 
             filename = file.path(path2figs, "study_arm_tsne.pdf"), 
             width = 15, height = 30)

## Process each subject separately
tsne_subject_res <- get_ordinations(ps, group = "Subject", method = "tsne")
saveRDS(tsne_subject_res, file = file.path(path2out, "subject_pca.rds"))

generate_pdf(ord_lst = pca_subject_res, 
             filename = file.path(path2figs, "subject_pca.pdf"), 
             width = 15, height = 25)



