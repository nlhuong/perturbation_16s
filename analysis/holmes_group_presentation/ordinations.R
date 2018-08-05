
#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## This is file contains code for performing ordination for 16S data.
##
## author: nlhuong90@gmail.com
## date: 07/11/2018

###############################################################################
## Setup and libraries
###############################################################################

rm(list = ls())
library("tidyverse")
library("phyloseq")
library("PMA")
library("adaptiveGPCA")


NOABX <- FALSE
AMPLICON <- FALSE
nPC <- 10

if(AMPLICON) {
  load("results/processed_physeq.rda")
  resfile <- "results/ordinate_16S.rda"
  #resfile <- "ordinate_16S_before_abx.rda"
  #load(resfile)
  
  ps <- ps_asinh
  if(NOABX) {
    ps <- subset_samples(ps, !(Group != "Abx") |  !grepl("MidAbx", Interval) | 
                           !grepl("PostAbx", Interval))
    norm_counts <- norm_counts[, sample_names(ps)]
  }
  
  taxtab <- data.frame(tax_table(ps)) %>%
    mutate(Seq_ID = rownames((.)))
  SMP <- sam_data(ps)
  norm_ihs <- t(asinh(norm_counts))
  
} else {
  load("results/metat_seed_filtered.rda")
  resfile <- "results/ordinate_metat.rda"
  taxtab <- data.frame(Seq_ID = rownames(covtab_seed), 
                       stringsAsFactors = FALSE) %>%
    separate(Seq_ID, into = paste0("SEED", 1:4), sep =";;", remove = FALSE)
  size_fac <- DESeq2::estimateSizeFactorsForMatrix(covtab_seed)
  covtab_seed_norm <- sweep(covtab_seed, MARGIN = 2, FUN = "/", size_fac)
  norm_ihs <- t(asinh(covtab_seed_norm))
  colnames(norm_ihs) <- taxtab$Seq_ID
  
  norm_ihs_mean <- data.frame(Subject = SMP$Subject, norm_ihs) %>%
    group_by(Subject) %>%
    summarise_all(mean) %>%
    as.data.frame() %>%
    column_to_rownames("Subject")
  norm_ihs_mean <- norm_ihs_mean[SMP$Subject, ]
  norm_ihs_centered <- as.matrix(norm_ihs) - as.matrix(norm_ihs_mean)
}


## PCA -------------------------------------------------------------------------
cat("Starting pca. \n")
ptm <- proc.time()
pca.ihs <- prcomp(norm_ihs, scale = FALSE)
time <- proc.time() - ptm
cat("Finished  pca in:\n")
time
# user  system elapsed 
# 8.568   0.012   8.580 


loadings <- pca.ihs$rotation[, 1:nPC] %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  rownames_to_column("Seq_ID")

if("taxtab" %in% ls()) {
  loadings <- loadings %>% left_join(taxtab)
}
scores <- pca.ihs$x[, 1:nPC] %>%
  as.data.frame() %>%
  rownames_to_column("Meas_ID") %>%
  left_join(SMP) %>%
  arrange(Group, Subject, Samp_Date)  %>%
  mutate(constant = "constant")

subject_scores <- scores %>%
  select(Subject, starts_with("PC")) %>%
  group_by(Subject) %>%
  summarise_all(mean)

scores_centered <- scores %>%
  left_join(subject_scores, by = c("Subject"), suffix = c("", ".subj")) 

for(PC in paste0("PC", 1:nPC)) {
  scores_centered[, PC] <- scores_centered[, PC] - 
    scores_centered[, paste0(PC, ".subj")]
}

cat("Starting centered pca. \n")
ptm <- proc.time()
pca.centered.ihs <- prcomp(norm_ihs_centered, scale = FALSE)
time <- proc.time() - ptm
cat("Finished  centered pca in:\n")
time
# user  system elapsed 
# 0.112   0.000   0.111

pca_centered_loadings <- data.frame(colnames(norm_ihs), 
                                    pca.centered.ihs$rotation[, 1:nPC])
colnames(pca_centered_loadings) <- 
  c("Seq_ID", paste0("PC_centered", seq_len(nPC)))

pca_centered_scores <- data.frame(rownames(norm_ihs), 
                                  pca.centered.ihs$x[, 1:nPC])
colnames(pca_centered_scores) <- 
  c("Meas_ID", paste0("PC_centered", seq_len(nPC)))

loadings <- loadings %>%
  left_join(pca_centered_loadings)

scores <- scores %>%
  left_join(pca_centered_scores) 


save(list = c("pca.ihs","pca.centered.ihs","scores", "loadings", "scores_centered"), 
     file = resfile)

## Sparse PCA ------------------------------------------------------------------
cat("Starting sparse pca. \n")
ptm <- proc.time()
sparse_pca <- PMA::SPC(
  scale(norm_ihs, center = TRUE, scale = FALSE),
  v = pca.ihs$rotation[, 1],
  K = 10, sumabsv = 15)

sparse_loadings <- data.frame(colnames(norm_ihs), sparse_pca$v)
colnames(sparse_loadings) <- 
  c("Seq_ID", paste0("sPC", seq_len(ncol(sparse_pca$v))))

sparse_scores <- data.frame(rownames(norm_ihs), sparse_pca$u)
colnames(sparse_scores) <- 
  c("Meas_ID", paste0("sPC",seq_len(ncol(sparse_pca$u))))

loadings <- loadings %>%
  left_join(sparse_loadings)

scores <- scores %>%
  left_join(sparse_scores) 

subject_scores <- scores %>%
  select(Subject,  starts_with("sPC")) %>%
  group_by(Subject) %>%
  summarise_all(mean)

scores_centered <- scores_centered %>%
  left_join(sparse_scores) %>%
  left_join(subject_scores, by = c("Subject"), suffix = c("", ".subj")) 

for(sPC in paste0("sPC", 1:5)) {
  scores_centered[, sPC] <- scores_centered[, sPC] - 
    scores_centered[, paste0(sPC, ".subj")]
}
save(list = c("pca.ihs", "pca.centered.ihs", "sparse_pca", "scores", "loadings", 
              "scores_centered"),
     file = resfile)
time <- proc.time() - ptm
cat("Finished sparse pca in:\n")
time
# user  system elapsed 
# 66.012   0.092  66.085 

## tSNE ------------------------------------------------------------------------

## Standard tSNE ----------------------------------
cat("Starting tsne. \n")
ptm <- proc.time()
set.seed(123)
rtsne.norm.ihs <- Rtsne::Rtsne(
  norm_ihs, dims = 2, initial_dims = 50, perplexity = 50,
  is_distance = FALSE, pca = TRUE, eta = 200, exaggeration_factor = 12)
time <- proc.time() - ptm
cat("Finished tsne in: \n")
time

## Centered tSNE ----------------------------------
nPC_for_tsne <- 50

centered_pca <- pca.ihs$x[, 1:nPC_for_tsne] %>%
  data.frame() %>%
  rownames_to_column("Meas_ID") %>%
  left_join(data.frame(SMP) %>% select(Meas_ID, Subject)) %>%
  group_by(Subject) %>%
  summarise_at(
    .vars = vars(contains("PC")),
    .funs = mean
  )

subjPCA <- data.frame(SMP) %>% 
  select(Meas_ID, Subject) %>%
  left_join(centered_pca) %>%
  select(-Subject) %>%
  column_to_rownames("Meas_ID")

centered_pca <- pca.ihs$x[, 1:nPC_for_tsne] - subjPCA[rownames(pca.ihs$x), ]

cat("Starting centered tsne. \n")
ptm <- proc.time()
set.seed(123)
rtsne.centered.ihs <- Rtsne::Rtsne(
  centered_pca, dims = 2, perplexity = 50, is_distance = FALSE, 
  pca = FALSE, eta = 200, exaggeration_factor = 12)
save(list = c("pca.ihs", "sparse_pca",
              "rtsne.norm.ihs", "rtsne.centered.ihs", "rtsne.bray.ihs", 
              "brayD.ihs", "pp", "out.agpca"),
     file = resfile)
time <- proc.time() - ptm
cat("Finished centered tsne in:\n")
time

## Bray tSNE ----------------------------------

cat("Computing bray distance")
ptm <- proc.time()
brayD.ihs <- vegan::vegdist(norm_ihs, method = "bray")
save(list = c("pca.ihs", "sparse_pca","rtsne.norm.ihs", "brayD.ihs"),
     file = resfile)
time <- proc.time() - ptm
cat("Finished bray dist in: \n")
time

cat("Starting bray tsne. \n")
ptm <- proc.time()
set.seed(123)
rtsne.bray.ihs <- Rtsne::Rtsne(
  brayD.ihs, is_distance = TRUE, dims = 2, perplexity = 50,
  pca = FALSE, eta = 200, exaggeration_factor = 12)
time <- proc.time() - ptm
cat("Finished bray tsne in:\n")
time

tsne_scores <- data.frame(
  Meas_ID = rownames(norm_ihs),
  tSNE1 = rtsne.norm.ihs$Y[, 1],
  tSNE2 = rtsne.norm.ihs$Y[, 2],
  tSNE1_bray = rtsne.bray.ihs$Y[, 1],
  tSNE2_bray = rtsne.bray.ihs$Y[, 2],
  tSNE1_centered = rtsne.centered.ihs$Y[, 1], 
  tSNE2_centered = rtsne.centered.ihs$Y[, 2]) 

scores <- scores %>%
  left_join(tsne_scores)

save(list = c("pca.ihs", "sparse_pca","rtsne.norm.ihs",
              "brayD.ihs", "rtsne.bray.ihs", "rtsne.centered.ihs",
              "scores", "loadings", "scores_centered"),
     file = resfile)


## agPCA -----------------------------------------------------------------------

if(AMPLICON) {
  cat("Starting agpca \n")
  ptm <- proc.time()
  pp <- processPhyloseq(ps)
  out.agpca <- adaptivegpca(pp$X, pp$Q, k = 5)
  
  agpca_scores <- data.frame(rownames(norm_ihs), out.agpca$U)
  colnames(agpca_scores) <- c("Meas_ID", paste0("agPC", 1:ncol(out.agpca$U)))
  
  agpca_loadings <- data.frame(colnames(norm_ihs), out.agpca$QV) 
  colnames(agpca_loadings) <- c("Seq_ID", paste0("agPC", 1:ncol(out.agpca$QV)))
  
  scores <- scores %>%
    left_join(agpca_scores) 
  
  loadings <- loadings %>%
    left_join(agpca_loadings) 
  
  save(list = c("pca.ihs", "sparse_pca",
                "rtsne.norm.ihs", "rtsne.centered.ihs", "rtsne.bray.ihs", 
                "brayD.ihs", 
                "pp", "out.agpca",
                "scores", "loadings", "scores_centered"),
       file = resfile)
  time <- proc.time() - ptm
  cat("Finished agpca in: \n")
  time
}




