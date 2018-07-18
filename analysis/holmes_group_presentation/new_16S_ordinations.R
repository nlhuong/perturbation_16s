
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

library("tidyverse")
library("phyloseq")
library("PMA")
library("adaptiveGPCA")

load("./processed_physeq.rda")
# load("ordinate_16S.rda")
# resfile <- "ordinate_16S.rda"
resfile <- "ordinate_16S_no_abx_cc.rda"
#resfile <- "ordinate_16S_before_abx.rda"


taxtab <- data.frame(tax_table(ps)) %>%
  mutate(Seq_ID = rownames((.)))

ps <- ps_asinh

### Can be deleted
ps <- subset_samples(ps,
  !(Group %in% c("Abx", "CC"))) # |  !grepl("MidAbx", Interval) |  !grepl("PostAbx", Interval))
#######################

norm_counts <- norm_counts[, sample_names(ps)]
norm_ihs <- t(asinh(norm_counts))


pca.ihs <- prcomp(norm_ihs, scale = FALSE)
save(list = c("pca.ihs"), file = resfile)

cat("Starting sparse pca. \n")
ptm <- proc.time()
sparse_pca <- PMA::SPC(
  scale(norm_ihs, center = TRUE, scale = FALSE),
  v = pca.ihs$rotation[, 1],
  K = 10, sumabsv = 10)
save(list = c("pca.ihs", "sparse_pca"), file = resfile)
time <- proc.time() - ptm
cat("Finished sparse pca in:\n")
time

cat("Starting tsne. \n")
ptm <- proc.time()
set.seed(123)
rtsne.norm.ihs <- Rtsne::Rtsne(
  norm_ihs, dims = 2, initial_dims = 50, perplexity = 50,
  is_distance = FALSE, pca = TRUE, eta = 200, exaggeration_factor = 12)
save(list = c("pca.ihs", "sparse_pca","rtsne.norm.ihs"),
     file = resfile)
time <- proc.time() - ptm
cat("Finished tsne in: \n")
time

cat("Computing bray distance")
ptm <- proc.time()
brayD.ihs <- phyloseq::distance(ps, method = "bray")
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
save(list = c("pca.ihs", "sparse_pca","rtsne.norm.ihs",
              "brayD.ihs", "rtsne.bray.ihs"),
     file = resfile)
time <- proc.time() - ptm
cat("Finished bray tsne in:\n")
time

cat("Starting agpca \n")
ptm <- proc.time()
pp <- processPhyloseq(ps)
out.agpca <- adaptivegpca(pp$X, pp$Q, k = 5)
save(list = c("pca.ihs", "sparse_pca","rtsne.norm.ihs",
              "brayD.ihs", "rtsne.bray.ihs", "pp", "out.agpca"),
     file = resfile)
time <- proc.time() - ptm
cat("Finished agpca in: \n")
time

centered_pca <- pca.ihs$x[, 1:50] %>%
  data.frame() %>%
  rownames_to_column("Meas_ID") %>%
  left_join(data.frame(sample_data(ps)) %>% 
              select(Meas_ID, Subject)) %>%
  group_by(Subject) %>%
  summarise_at(
    .vars = vars(contains("PC")),
    .funs = mean
  )

subjPCA <- data.frame(sample_data(ps)) %>% 
  select(Meas_ID, Subject) %>%
  left_join(centered_pca) %>%
  select(-Subject) %>%
  column_to_rownames("Meas_ID")

centered_pca <- pca.ihs$x[, 1:50] - subjPCA[rownames(pca.ihs$x), ]

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

################################################################################

taxtab <- data.frame(tax_table(ps)) %>%
  mutate(Seq_ID = rownames((.)))

nPC <- 5
loadings <- pca.ihs$rotation[, 1:nPC] %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  rownames_to_column("Seq_ID") %>%
  left_join(taxtab)

scores <- pca.ihs$x[, 1:nPC] %>%
  as.data.frame() %>%
  rownames_to_column("Meas_ID") %>%
  left_join(SMP) %>%
  arrange(Group, Subject, Samp_Date)  %>%
  mutate(constant = "constant")

subject_scores <- scores %>%
  select(Subject, PC1:PC5) %>%
  group_by(Subject) %>%
  summarise_all(mean)

scores_centered <- scores %>%
  left_join(subject_scores, by = c("Subject"), suffix = c("", ".subj")) %>%
  mutate(
    PC1 = PC1 - PC1.subj,
    PC2 = PC2 - PC2.subj,
    PC3 = PC3 - PC3.subj,
    PC4 = PC4 - PC4.subj,
    PC5 = PC5 - PC5.subj
  )

sparse_loadings <- data.frame(taxa_names(ps), sparse_pca$v)
colnames(sparse_loadings) <- c("Seq_ID", paste0("sPC", 1:ncol(sparse_pca$v)))

sparse_scores <- data.frame(sample_names(ps), sparse_pca$u)
colnames(sparse_scores) <- c("Meas_ID", paste0("sPC", 1:ncol(sparse_pca$u)))

loadings <- loadings %>%
  left_join(sparse_loadings)

scores <- scores %>%
  left_join(sparse_scores) 


agpca_scores <- data.frame(sample_names(ps), out.agpca$U)
colnames(agpca_scores) <- c("Meas_ID", paste0("agPC", 1:ncol(out.agpca$U)))

agpca_loadings <- data.frame(taxa_names(ps), out.agpca$QV) 
colnames(agpca_loadings) <- c("Seq_ID", paste0("agPC", 1:ncol(out.agpca$QV)))

scores <- scores %>%
  left_join(agpca_scores) 

loadings <- loadings %>%
  left_join(agpca_loadings) 

tsne_scores <- data.frame(
  Meas_ID = sample_names(ps),
  tSNE1 = rtsne.norm.ihs$Y[, 1],
  tSNE2 = rtsne.norm.ihs$Y[, 2],
  tSNE1_bray = rtsne.bray.ihs$Y[, 1],
  tSNE2_bray = rtsne.bray.ihs$Y[, 2],
  tSNE1_centered = rtsne.centered.ihs$Y[, 1], 
  tSNE2_centered = rtsne.centered.ihs$Y[, 2]) 

scores <- scores %>%
  left_join(tsne_scores)

save(list = c("pca.ihs", "sparse_pca",
              "rtsne.norm.ihs", "rtsne.centered.ihs", 
              "brayD.ihs", "rtsne.bray.ihs", 
              "pp", "out.agpca",
              "scores", "loadings", "scores_centered"), 
     file = resfile)






