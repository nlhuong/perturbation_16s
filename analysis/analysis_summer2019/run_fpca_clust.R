library(phyloseq)
library(fdapace)
library(tidyverse)
source("fpca_funs.R")
NCORES <- min(16, parallel::detectCores())
RUN_ABX <- FALSE
RUN_DIET <- TRUE


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
  select(-Seq) %>%
  mutate(
    OrgName.RDP = paste(Genus, Species),
    OrgName.Silva = paste(GenusSilva, SpeciesSilva),
    OrgName = ifelse(grepl("NA", OrgName.RDP) & !grepl("NA", OrgName.RDP),
                     OrgName.Silva, OrgName.RDP),
    OrgName = ifelse(OrgName == "NA NA", "NA", OrgName),
    OrgName = paste0(Seq_ID, ": ", OrgName)) %>%
  select(1:8, OrgName) 
head(TAXTAB)

## Process data
dds <- phyloseq_to_deseq2(psSubj, design = ~ 1)
dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
norm_counts <- DESeq2::counts(dds, normalized = TRUE)
psSubj.norm <- psSubj
otu_table(psSubj.norm) <- otu_table(norm_counts, taxa_are_rows = TRUE)
psSubj.norm

if(RUN_ABX) {

    abxFac <- c("PreAbx", "MidAbx", "PostAbx")
    abx <- subset_samples(psSubj.norm, Abx_RelDay >= -50 & Abx_RelDay <= 60)
    abx <- subset_taxa(abx, taxa_sums(abx) > 0)
    print(abx)

    (replicated <- data.frame(sample_data(abx)) %>% 
        group_by(Subject, Group, Interval, Abx_RelDay) %>%
        mutate(n = n()) %>%
         ungroup() %>% filter(n > 1))
    print(replicated)
    # Remove replicated samples
    print(sample_sums(abx)[replicated$Meas_ID])

    #  M2242  M2243  M7869  M7886 
    # 109520  88834  74697  74873 

    # we retain a replicate with higher sample depth:

    abx <- subset_samples(abx, !Meas_ID %in% c("M2243", "M7869"))
    print(abx)

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

    cat("fitting", length(unique(keep_rsv$Seq_ID)), "sequences\n") # = 920 

    abx.rsv.subset <- abx.rsv %>%
        filter(Seq_ID %in% unique(keep_rsv$Seq_ID))

    save(list = c("keep_rsv", "abx.rsv.subset"), 
         file = "output/fpca_clust_res_abx_norm.rda")

    abx_fclust <- fpca_wrapper(
        abx.rsv.subset,
        time_column = "Abx_RelDay",
        value_column = "abundance",
        replicate_column = "Subject",
        feat_column = "Seq_ID",
        cluster = TRUE,
        clust_min_num_replicate = 15,
        fpca_optns = NULL, fclust_optns = NULL,
        parallel = TRUE, ncores = NCORES)


    save(list = c("keep_rsv", "abx.rsv.subset","abx_fclust"),
        file = "output/fpca_clust_res_abx_norm.rda")

    abx_fclust_mu <- get_fpca_means(abx_fclust)
    abx_fclust_subjFit <- get_fpca_fits(abx_fclust)

    save(list = c("keep_rsv", "abx.rsv.subset", 
                 "abx_fclust", "abx_fclust_subjFit", "abx_fclust_mu"),
        file = "output/fpca_clust_res_abx_norm.rda")
    
    set.seed(123456)
    min_num_subj <- 5
    
    abx2plot.fclust <- abx_fclust_subjFit %>%
      left_join(TAXTAB, by = c("Feature_ID" = "Seq_ID")) %>%
      group_by(Feature_ID) %>%
      mutate(
        n_cluster1 = sum(Cluster == 1),
        n_cluster2 = sum(Cluster == 2),
        majorityCluster = ifelse(n_cluster1 > n_cluster2, 1, 2)
      ) %>% 
      ungroup() %>%
      mutate(
        Subject_Cluster = ifelse(Cluster == majorityCluster, "Majority", "Minority")) 
    
    abx_fclust_majority <- abx2plot.fclust %>%
      filter(Subject_Cluster == "Majority") %>%
      left_join(keep_rsv, by = c("Feature_ID" = "Seq_ID")) %>%
      filter(Freq >= min_num_subj) %>%
      select(-Freq) %>%
      group_by(time, Feature_ID) %>%
      summarise(value = mean(value)) %>%
      left_join(TAXTAB, by = c("Feature_ID" = "Seq_ID")) %>%
      arrange( Feature_ID, time)
    
    length(unique(abx_fclust_majority$Feature_ID))
    #[1] 573
    
    seqClusters.abx <- fit_fpca(
      abx_fclust_majority,
      "time", "value", "Feature_ID",
      cluster = TRUE, K = 8)
    
    seqClusters_fpca.abx <- fitted_values_fpca(seqClusters.abx) %>%
      mutate(Seq_ID = Replicate_ID) 
    
    save(list = c("keep_rsv", "abx.rsv.subset", 
                  "abx_fclust", "abx_fclust_subjFit", "abx_fclust_mu",
                  "abx_fclust_majority", "seqClusters.abx", "seqClusters_fpca.abx"),
         file = "output/fpca_clust_res_abx_norm.rda")
}

##-------------------------------------------------------------------------------------------

if(RUN_DIET) {

    dietFac <- c("PreDiet", "MidDiet", "PostDiet")
    diet <- subset_samples(psSubj.norm, Diet_RelDay >= -30 & Diet_RelDay <= 45)
    diet <- subset_taxa(diet, taxa_sums(diet) > 0)
    print(diet)


    (replicated <- data.frame(sample_data(diet)) %>% 
        group_by(Subject, Group, Diet_Interval, Diet_RelDay) %>%
        mutate(n = n()) %>% ungroup() %>% 
        filter(n > 1) %>% select(Meas_ID, Subject, Diet_RelDay, Diet_Interval))
    print(replicated)

    # Remove replicated samples
    print(sample_sums(diet)[replicated$Meas_ID])


    # we retain a replicate with higher sample depth:
    diet <- subset_samples(diet, !Meas_ID %in% c("M2129", "M7869","M8253"))
    print(diet)

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

    cat("fitting", length(unique(keep_rsv$Seq_ID)), "sequences\n") 

    diet.rsv.subset <- diet.rsv %>%
        filter(Seq_ID %in% unique(keep_rsv$Seq_ID))

    save(list = c("keep_rsv", "diet.rsv.subset"), 
         file = "output/fpca_clust_res_diet_norm.rda")


    diet_fclust <- fpca_wrapper(
        diet.rsv.subset,
        time_column = "Diet_RelDay",
        value_column = "abundance",
        replicate_column = "Subject",
         feat_column = "Seq_ID",
        cluster = TRUE,
        clust_min_num_replicate = 15,
        fpca_optns = NULL, fclust_optns = NULL,
        parallel = TRUE, ncores = NCORES)


    save(list = c("keep_rsv", "diet.rsv.subset","diet_fclust"),
        file = "output/fpca_clust_res_diet_norm.rda")

    diet_fclust_mu <- get_fpca_means(diet_fclust)
    diet_fclust_subjFit <- get_fpca_fits(diet_fclust)

    save(list = c("keep_rsv", "diet.rsv.subset", 
              "diet_fclust", "diet_fclust_subjFit", "diet_fclust_mu"),
        file = "output/fpca_clust_res_diet_norm.rda")
    
    set.seed(123456)
    min_num_subj <- 5

    diet2plot.fclust <- diet_fclust_subjFit %>%
      left_join(TAXTAB, by = c("Feature_ID" = "Seq_ID")) %>%
      group_by(Feature_ID) %>%
      mutate(
        n_cluster1 = sum(Cluster == 1),
        n_cluster2 = sum(Cluster == 2),
        majorityCluster = ifelse(n_cluster1 > n_cluster2, 1, 2)
      ) %>% 
      ungroup() %>%
      mutate(
        Subject_Cluster = ifelse(Cluster == majorityCluster, "Majority", "Minority")) 
    
    diet_fclust_majority <- diet2plot.fclust %>%
      filter(Subject_Cluster == "Majority") %>%
      left_join(keep_rsv, by = c("Feature_ID" = "Seq_ID")) %>%
      filter(Freq >= min_num_subj) %>%
      select(-Freq) %>%
      group_by(time, Feature_ID) %>%
      summarise(value = mean(value)) %>%
      left_join(TAXTAB, by = c("Feature_ID" = "Seq_ID")) %>%
      arrange( Feature_ID, time)
    
    length(unique(diet_fclust_majority$Feature_ID))
    #[1] 500
    
    set.seed(123456)
    min_num_subj <- 5
    seqClusters.diet <- fit_fpca(
      diet_fclust_majority,
      "time", "value", "Feature_ID",
      cluster = TRUE, K = 6)
    
    seqClusters_fpca.diet <- fitted_values_fpca(seqClusters.diet) %>%
      mutate(Seq_ID = Replicate_ID) 
    
    save(list = c("keep_rsv", "diet.rsv.subset", 
                  "diet_fclust", "diet_fclust_subjFit", "diet_fclust_mu",
                  "diet_fclust_majority", "seqClusters.diet", "seqClusters_fpca.diet"),
         file = "output/fpca_clust_res_diet_norm.rda")
    
}
