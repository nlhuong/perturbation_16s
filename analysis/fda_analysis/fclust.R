setwd("/scratch/users/lanhuong/Projects/perturbation_16s/analysis/fda_analysis")
NCORES <- 16

library("DESeq2")
library("phyloseq")
library("fdapace")
library("viridis")
library("RColorBrewer")
library("ggplot2")
library("tidyverse")
source('../utils.R')
source('./fda_pca_funs.R')

datadir <- "/scratch/PI/sph/resilience/"
load("results/ps_objects.rda")
load("results/abx_long_subset.rda")

# duplicated_meas <- data.frame(abx@sam_data) %>%
#   select(Meas_ID, Subject, Abx_RelDay, Samp_Date, Samp_Time) %>%
#   group_by(Subject, Abx_RelDay) %>%
#   filter(n() > 1)
# duplicated_meas
# 
# abx_SMP <- data.frame(abx@sam_data) %>%
#   mutate(
#     Abx_RelDay = ifelse(Meas_ID == "M2243", 68.5, Abx_RelDay),
#     Abx_RelDay = ifelse(Meas_ID == "M7886", -48.5, Abx_RelDay))
# rownames(abx_SMP) <- abx_SMP$Meas_ID 
# sample_data(abx) <- sample_data(abx_SMP)
# 
# abx.rsv <- data.frame(asinh(as(otu_table(abx), "matrix"))) %>%
#   rownames_to_column("Seq_ID") %>%
#   gather(Meas_ID, abundance, -Seq_ID) %>%
#   left_join(abx_SMP) %>%
#   left_join(taxtab) %>%
#   mutate(
#     Abx_Interval = factor(
#       Abx_Interval, levels = c("PreAbx", "MidAbx", "PostAbx"))) %>%
#   filter(Abx_RelDay >= -40, Abx_RelDay <= 80) %>%   # limit the scope due to effect of other perturbations
#   arrange(Seq_ID, Subject, Abx_RelDay)
# 
# thresh <- 0; num_nonzero <- 10; min_no_subj <- 3
# keep_rsv <- abx.rsv %>%
#   filter(abundance > thresh) %>%
#   group_by(Subject, Seq_ID) %>%
#   summarise(Freq = n()) %>%   # No. of samples for each subject with non-zero count of a given rsv
#   filter(Freq >= num_nonzero) %>% 
#   group_by(Seq_ID) %>%
#   summarise(Freq = n()) %>% # For each RSV, # of subjects w/ num_nonzero samples with counts > thresh
#   filter(Freq >= min_no_subj) %>% 
#   arrange(desc(Freq))
# 
# abx.rsv.subset <- abx.rsv %>% 
#   filter(Seq_ID %in% unique(keep_rsv$Seq_ID))

length(unique(keep_rsv$Seq_ID)) # = 861


# abx_fpca <- fda_pca(
#   abx.rsv.subset,
#   time_column = "Abx_RelDay",
#   value_column = "abundance",
#   replicate_column = "Subject",
#   feat_column = "Seq_ID",
#   cluster = FALSE,
#   fpca_optns = NULL, fclust_optns = NULL,
#   parallel = TRUE, ncores = NCORES)
# abx_fpca_mu <- get_fpca_means(abx_fpca)
# abx_fpca_subjFit <- get_fpca_fits(abx_fpca)
# save(list = c("abx_fpca", "abx_fpca_mu", "abx_fpca_subjFit"),
#      file = "results/abx_rsv_fpca.rda")


abx_fclust <- fda_pca(
  abx.rsv.subset,
  time_column = "Abx_RelDay",
  value_column = "abundance",
  replicate_column = "Subject",
  feat_column = "Seq_ID",
  cluster = TRUE,
  fpca_optns = NULL, fclust_optns = NULL,
  parallel = FALSE, ncores = NCORES)
save(list = c("abx_fclust"),
     file = "results/abx_rsv_fclust.rda")

abx_fclust_mu <- get_fpca_means(abx_fclust)
abx_fclust_subjFit <- get_fpca_fits(abx_fclust)

save(list = c("abx_fclust", "abx_fclust_subjFit", "abx_fclust_mu"),
     file = "results/abx_rsv_fclust.rda")
