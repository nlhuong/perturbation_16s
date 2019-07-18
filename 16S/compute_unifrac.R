library(phyloseq)


psSubj <- readRDS("../data/16S/phyloseq/perturb_physeq_participants_decontam_15Jul19.rds")
psSubj

load("../analysis/analysis_summer2019/output/pairwise_dist_subj_16S.rda")

study_arms <- as.character(unique(psSubj@sam_data$Group))

uniFracD.lst <- list()
study_arms <- setdiff(study_arms, names(uniFracD.lst))
for(gr in study_arms) {
  ps.gr <- phyloseq::subset_samples(psSubj, Group == gr)
  print(ps.gr)
  uniFracD.lst[[gr]] <- phyloseq::distance(
    ps.gr, method = "unifrac", parallel = TRUE, 
    fast = TRUE, normalized = TRUE, weighted = FALSE)
  save(list = c("brayD.ihs", "jaccardD", "uniFracD.lst"), 
       file = "../analysis/analysis_summer2019/output/pairwise_dist_subj_16S.rda")
}
 



