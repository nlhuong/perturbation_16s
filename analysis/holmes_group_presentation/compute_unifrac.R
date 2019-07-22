library("phyloseq")

load( "results/ordinate_16S.rda")
load("results/processed_physeq.rda")

subject_list <- unique(sample_data(ps)$Subject)
cat("Number of subjects: ", length(subject_list), ".\n")
unifrac_list <- list()
i <- 0
for(subj in subject_list) {
    i <- i +1
    cat("Computing UniFrac dist for Subject: ", subj, ".\n")
    subj_ps <- subset_samples(ps, Subject == subj)
    subj_ps <- subset_taxa(subj_ps, taxa_sums(subj_ps) > 0)
    doParallel::registerDoParallel(cores = 16)
    subj_unifracD <- phyloseq::UniFrac(subj_ps, parallel=FALSE, fast=TRUE)
    unifrac_list[[subj]] <- reshape2::melt(
        as.matrix(subj_unifracD), varnames = c("S1", "S2"),
        mutate(S1 = as.character(S1),
               S2 = as.character(S2))
    )
}
unifrac_df <- plyr::ldply(unifrac_list, .id = "Subject")
colnames(unifrac_df) <- c("Subject", "S1", "S2", "unifrac")
save(list = c("unifrac_list", "unifrac_df"), file = "results/unifrac_dist.rda")


