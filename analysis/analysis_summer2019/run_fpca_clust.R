library(phyloseq)
library(fdapace)
library(tidyverse)
source("fpca_funs.R")

RUN_ABX <- TRUE
RUN_DIET <- TRUE
RUN_CC <- TRUE


# File generated in /perturbation_16s/analysis/analysis_summer2019/generate_phyloseq.rmd
psSubj <- readRDS(
  "../../data/16S/phyloseq/perturb_physeq_participants_decontam_15Jul19.rds")
psSubj
# otu_table()   OTU Table:         [ 2425 taxa and 4402 samples ]
# sample_data() Sample Data:       [ 4402 samples by 40 sample variables ]
SMP <- data.frame(sample_data(psSubj))
SUBJ <- SMP %>% select(Subject, Age:BirthYear) %>% distinct()
TAXTAB <- data.frame(as(tax_table(psSubj), "matrix")) %>%
  rownames_to_column("Seq_ID") %>%
  select(1:7) 

if(RUN_ABX) {

    abxFac <- c("PreAbx", "MidAbx", "PostAbx")
    abx <- subset_samples(psSubj, Abx_RelDay >= -50, Abx_RelDay <= 60)
    abx <- subset_taxa(abx, taxa_sums(abx) > 0)
    abx

    (replicated <- data.frame(sample_data(abx)) %>% 
        group_by(Subject, Group, Interval, Abx_RelDay) %>%
        mutate(n = n()) %>%
         ungroup() %>% filter(n > 1))

    # Remove replicated samples
    sample_sums(abx)[replicated$Meas_ID]

    #  M2242  M2243  M7869  M7886 
    # 109520  88834  74697  74873 

    # we retain a replicate with higher sample depth:

    abx <- subset_samples(abx, !Meas_ID %in% c("M2243", "M7869"))
    abx

    abx.rsv <- data.frame(asinh(as(otu_table(abx), "matrix"))) %>%
        rownames_to_column("Seq_ID") %>%
        gather(Meas_ID, abundance, -Seq_ID) %>%
        left_join(SMP) %>%
        left_join(TAXTAB) %>%
        mutate(Abx_Interval = factor(Abx_Interval, levels = abxFac)) %>%
        arrange(Seq_ID, Subject, Abx_RelDay)


    thresh <- 0; num_nonzero <- 10; 
    min_no_subj <- 3

    keep_rsv <- abx.rsv %>%
        filter(abundance > thresh) %>%
        group_by(Subject, Seq_ID) %>%
        summarise(Freq = n()) %>%   # No. of samples per subject w/ positive rsv count 
        filter(Freq >= num_nonzero) %>%
        group_by(Seq_ID) %>%
        summarise(Freq = n()) %>% # Number of subjects w/ num_nonzero samples with counts > thresh
        filter(Freq >= min_no_subj) %>%
        arrange(desc(Freq))

    length(unique(keep_rsv$Seq_ID)) # = 920

    abx.rsv.subset <- abx.rsv %>%
        filter(Seq_ID %in% unique(keep_rsv$Seq_ID))

    save(list = c("keep_rsv", "abx.rsv.subset"), 
         file = "output/fpca_clust_res_abx.rda")

    abx_fclust <- fpca_wrapper(
        abx.rsv.subset,
        time_column = "Abx_RelDay",
        value_column = "abundance",
        replicate_column = "Subject",
        feat_column = "Seq_ID",
        cluster = TRUE,
        clust_min_num_replicate = 15,
        fpca_optns = NULL, fclust_optns = NULL,
        parallel = TRUE, ncores = 16)


    save(list = c("keep_rsv", "abx.rsv.subset","abx_fclust"),
        file = "output/fpca_clust_res_abx.rda")

    abx_fclust_mu <- get_fpca_means(abx_fclust)
    abx_fclust_subjFit <- get_fpca_fits(abx_fclust)

    save(list = c("keep_rsv", "abx.rsv.subset", 
                 "abx_fclust", "abx_fclust_subjFit", "abx_fclust_mu"),
        file = "output/fpca_clust_res_abx.rda")
}

##-------------------------------------------------------------------------------------------

if(RUN_DIET) {

    dietFac <- c("PreDiet", "MidDiet", "PostDiet")
    diet <- subset_samples(psSubj, Diet_RelDay >= -30, Diet_RelDay <= 30)
    diet <- subset_taxa(diet, taxa_sums(diet) > 0)
    diet


    (replicated <- data.frame(sample_data(diet)) %>% 
        group_by(Subject, Group, Diet_Interval, Diet_RelDay) %>%
        mutate(n = n()) %>% ungroup() %>% 
        filter(n > 1) %>% select(Meas_ID, Subject, Diet_RelDay, Diet_Interval))

    # Remove replicated samples
    sample_sums(diet)[replicated$Meas_ID]


    # we retain a replicate with higher sample depth:
    diet <- subset_samples(diet, !Meas_ID %in% c("M2129", "M7869","M8253"))
    diet

    diet.rsv <- data.frame(asinh(as(otu_table(diet), "matrix"))) %>%
        rownames_to_column("Seq_ID") %>%
        gather(Meas_ID, abundance, -Seq_ID) %>%
        left_join(SMP) %>%
        left_join(TAXTAB) %>%
        mutate(Diet_Interval = factor(Diet_Interval, levels = dietFac)) %>%
        arrange(Seq_ID, Subject, Diet_RelDay)

    thresh <- 0; num_nonzero <- 10; 
    min_no_subj <- 3

    keep_rsv <- diet.rsv %>%
        filter(abundance > thresh) %>%
        group_by(Subject, Seq_ID) %>%
        summarise(Freq = n()) %>%   # No. of samples per subject w/ positive rsv count 
        filter(Freq >= num_nonzero) %>%
        group_by(Seq_ID) %>%
        summarise(Freq = n()) %>% # Number of subjects w/ num_nonzero samples with counts > thresh
        filter(Freq >= min_no_subj) %>%
        arrange(desc(Freq))

    length(unique(keep_rsv$Seq_ID)) 

    diet.rsv.subset <- diet.rsv %>%
        filter(Seq_ID %in% unique(keep_rsv$Seq_ID))

    save(list = c("keep_rsv", "diet.rsv.subset"), 
         file = "output/fpca_clust_res_diet.rda")


    diet_fclust <- fpca_wrapper(
        diet.rsv.subset,
        time_column = "Diet_RelDay",
        value_column = "abundance",
        replicate_column = "Subject",
         feat_column = "Seq_ID",
        cluster = TRUE,
        clust_min_num_replicate = 15,
        fpca_optns = NULL, fclust_optns = NULL,
        parallel = TRUE, ncores = 16)


    save(list = c("keep_rsv", "diet.rsv.subset","diet_fclust"),
        file = "output/fpca_clust_res_diet.rda")

    diet_fclust_mu <- get_fpca_means(diet_fclust)
    diet_fclust_subjFit <- get_fpca_fits(diet_fclust)

    save(list = c("keep_rsv", "diet.rsv.subset", 
              "diet_fclust", "diet_fclust_subjFit", "diet_fclust_mu"),
        file = "output/fpca_clust_res_diet.rda")
}
