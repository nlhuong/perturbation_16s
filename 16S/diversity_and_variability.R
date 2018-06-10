library(phyloseq)
library(tidyverse)
rm(list = ls())
theme_set(theme_bw())


path_to_data <- "perturbation_16s/data/processed/"
physeq_file <- file.path(path_to_data, "perturb_physeq_filtered_22Jan18.rds")
load(file.path(path_to_data, "perturb_ordinate.rda"))
ps0 <- readRDS(physeq_file)
ps0
theme_set(theme_bw())

remove_subjects <- c("AAA", "AAB", "AAN", "DAC")
ps <- subset_samples(ps0, !Subject %in% remove_subjects)
ps

# Filter taxa present in more than 1 subject
minTaxaSubjPrev <- 2
# Compute prevalence across subjects
num_subj_present <- function(physeq, subject_indicator = "Subject") {
  smp <- data.frame(physeq@sam_data, stringsAsFactors = FALSE)
  seqtab <- as(physeq@otu_table, "matrix")
  sapply(1:ntaxa(physeq), function(i) {
    subj.present <- smp[[subject_indicator]][seqtab[i, ] > 0]
    length(unique(subj.present))})
}
subjects_prevalence <- num_subj_present(ps)
ps <- subset_taxa(ps, subjects_prevalence >= minTaxaSubjPrev)
ps


# SMP <- data.frame(ps@sam_data, stringsAsFactors = FALSE) %>%
#   mutate(
#     SampleDepth = sample_sums(ps),
#     Interval = paste0(Diet_Interval, "_", CC_Interval, "_", Abx_Interval),
#     Interval = gsub("_NA", "", Interval),
#     Interval = gsub("NA_", "", Interval),
#     Study_Arm = ifelse(Group == "CC_Diet", "Diet_CC", Group),
#     Study_Arm = ifelse(Study_Arm == "CC_Diet_Abx", "Diet_Abx_CC", Study_Arm),
#     Study_Arm = factor(Study_Arm, levels = c("NoIntv", "Abx", "Diet_Abx", "Diet_CC", "Diet_Abx_CC"))
# ) 
# 
# SUBJ <-  SMP %>%
#   select(Subject, Group, First_Sample_Date, CC_Date, Diet_StartDate, 
#          Age, Gender, Height, Weight, BMI,
#          Waist_Circumference, Resting_Pulse, Blood_Pressure,
#          Test_Cholesterol_Glucose, Birth_Mode, BirthYear) %>%
#   distinct()


################################################################################

brayD <- as(brayD, "matrix") 
colnames(brayD) <- rownames(brayD) <- factor(rownames(brayD), ordered = TRUE)

DF <- ps@sam_data %>%
  as.data.frame()

brayDF <- brayD %>%
  reshape2::melt(value.name = "bray_dist", varnames = c("i", "j")) %>%
  as.data.frame() %>%
  mutate(
    i = factor(i, level = rownames(brayD), ordered = TRUE),
    j = factor(j, level = rownames(brayD), ordered = TRUE)
  ) %>%
  left_join(DF %>%select(1:20), by = c("i" = "Meas_ID")) %>%
  left_join(DF %>% select(1:20), by = c("j" = "Meas_ID")) %>%
  filter(Subject.x == Subject.y, i != j) 

brayDF1 <- brayDF %>%
  mutate(
    day_diff = as.numeric(Samp_Date.y - Samp_Date.x),
    abs_day_diff = abs(day_diff),
    thresh_day_diff = pmin(abs_day_diff, 5),
    bray_per_day = bray_dist/thresh_day_diff
  ) 

idx <- brayDF1 %>%
  group_by(i) %>%
  summarise(
    min_diff_day_before = ifelse(sum(day_diff <= 0) >= 1,
                                 max(day_diff[day_diff <= 0]),
                                 -Inf),
    min_diff_day_after =  ifelse(sum(day_diff > 0) >= 1,
                                 min(day_diff[day_diff > 0]),
                                 +Inf)
  ) %>%
  filter(abs(min_diff_day_after) > 1 | abs(min_diff_day_before) > 1)

brayDF.g <- brayDF1 %>%
  group_by(i) %>%
  top_n(5, wt = -day_diff)

brayDF.m <- brayDF.g %>%
  group_by(i, Group.x, Subject.x, Timeline.x, DayFromStart.x, 
           Samp_Type.x, Abx_Interval.x ) %>% 
  summarise(
    mean_bray_per_day = mean(bray_per_day),
    mean_bray = mean(bray_dist),
    median_bray_per_day = median(bray_per_day),
    median_bray = median(bray_dist)
  )

ggplot(brayDF.g %>% filter(Group.x == "Diet_Abx"),
       aes(x = DayFromStart.x, y = bray_per_day, color = Timeline.x)
) +
  geom_point(
    aes(y = bray_per_day),
    alpha = 0.3, size = 1
  ) +
  geom_line(
    data = brayDF.m %>% filter(Group.x == "Diet_Abx"),
    aes(y = mean_bray_per_day, group = Subject.x), 
    size = 1
  ) +
  geom_point(
    data = brayDF.m %>% filter(Group.x == "Diet_Abx"),
    aes(y = mean_bray_per_day),
    size = 1
  ) +
  facet_wrap(~ Subject.x, ncol = 2, scales = "free") +
  scale_color_manual(values = c("#F8766D", "#00BFC4", "grey60")) +
  #ylim(NA, 0.8) +
  theme(text = element_text(size = 20))



###############################################################################

psgenus <- tax_glom(ps0, "Genus", NArm = TRUE)
psgenus <- readRDS(file.path(path_to_data, "perturb_physeq_genus.rds"))
psgenus
saveRDS(psgenus, 
        file = file.path(path_to_data, "perturb_physeq_genus.rds"))
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 197 taxa and 2965 samples ]
# sample_data() Sample Data:       [ 2965 samples by 270 sample variables ]
# tax_table()   Taxonomy Table:    [ 197 taxa by 12 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 197 tips and 196 internal nodes ]



H <- as(psgenus@otu_table, "matrix") + 0.1
H <- sweep(H, 2, colSums(H), "/")
H <- -colSums(H*log(H))

Richness <- colSums(psgenus@otu_table > 0)

SMP <- data.frame(psgenus@sam_data) %>%
  mutate(
    Richness = Richness,
    Shannon = H
  )

plt <- ggplot(
  SMP %>% 
    filter(Group == "Diet_Abx", 
           Subject %in% c("EAY", "EAQ", "EAP", "EAR", "EAS", "EAT" )),
  aes(x = DayFromStart, y = Richness, group = Subject, color = Timeline)
  ) +
  geom_line(alpha = 0.8, size = 0.7) +
  geom_point() +
  geom_line(aes(y = Shannon*10), alpha = 0.5, size = 0.7) +
  geom_point(aes(y=  Shannon*10), pch = 17) +
  facet_wrap(~ Subject, ncol = 1, scales = "free_x") +
  scale_y_continuous(name = "Richness [Distinct Genus Count]", 
                     sec.axis = sec_axis(~./10, name = "Shannon Index")) +
  scale_color_manual(values = c("#F8766D", "#00BFC4", "grey60"))
  

ggsave("perturbation_16s/figs/diversity_diet_abx.pdf", 
       plt,
       width = 15, height = 10)

plt <- ggplot(
  SMP %>% 
    filter(Group == "NoIntv", 
           Subject %in% c("CAA", "CAC", "CAL", "CAD", "CAM","CAP")), 
  #c("NA_QC", "DAC")),
  aes(x = DayFromStart, y = Richness, group = Subject, color = Abx_Interval)
  ) +
  geom_line(alpha = 0.8, size = 0.7) +
  geom_point() +
  geom_line(aes(y = Shannon*10), alpha = 0.5, size = 0.7) +
  geom_point(aes(y=  Shannon*10), pch = 17) +
  facet_wrap(~ Subject, ncol = 1, scales = "free_x") +
  scale_y_continuous(name = "Richness [Distinct Genus Count]", 
                     sec.axis = sec_axis(~./10, name = "Shannon Index")) +
  scale_color_manual(name = "Timeline", values = c("grey60", "#F8766D"))

ggsave("perturbation_16s/figs/diversity_nointv.pdf", 
       plt,
       width = 15, height = 10)

plt <- ggplot(SMP %>% filter(Group == "Abx")) +
  geom_line(aes(x = DayFromStart, y = Richness, group = Subject, 
                color = Timeline), size = 1) +
  facet_grid(Subject ~ Group) +
  scale_color_manual(values = c("#F8766D", "grey60"))

ggsave("perturbation_16s/figs/richness_abx.pdf", 
       plt,
       width = 10, height = 12)
