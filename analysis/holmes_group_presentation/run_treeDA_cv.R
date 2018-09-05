rm(list = ls())
library("dplyr")
library("phyloseq")
library("treeDA")
library("adaptiveGPCA")
library("DESeq2")

datadir <- "/scratch/PI/sph/resilience/"
curdir <- getwd()
resfile <- "results/treeDA_res16S.rda"
source('utils.R')

processed_datafile <- "results/processed_physeq.rda"
load(processed_datafile)

###############################################################################
## Process Data
###############################################################################
# Normalize by size factors
dds <- phyloseq_to_deseq2(ps, design = ~ 1)
dds <- estimateSizeFactors(dds, type = "poscounts")
norm_counts <- counts(dds, normalized = TRUE)
ps_norm <- ps
otu_table(ps_norm) <- otu_table(norm_counts, taxa_are_rows = TRUE)

## Diet Samples ----------------------------------------------------------------
# Subest to only diet pre/post diet but before other perturbations
# Subest to only diet pre/post diet but before other perturbations
diet <- subset_samples(ps_norm, grepl("Diet", Group))
diet <- subset_samples(
  diet, !grepl("MidAbx", Abx_Interval) & !grepl("PostAbx", Abx_Interval) &
    !grepl("PostCC", CC_Interval))
diet

# Filter taxa that occurs at minTaxaSum in at least minNoSubj subjects
minTaxaSum <- 5; minTaxaPrev <- 3; minNoSubj <- 3
ASVsum <- data.frame(t(as(otu_table(diet), "matrix"))) %>%
  mutate(Subject = diet@sam_data$Subject) %>%
  group_by(Subject) %>%
  summarise_all(sum)
ASVprev <- data.frame(t(as(otu_table(diet), "matrix")) > 0) %>%
  mutate(Subject = diet@sam_data$Subject) %>%
  group_by(Subject) %>%
  summarise_all(sum)

keepASV <- data.frame(
  Seq_ID = colnames(ASVsum),
  abundant = colSums(ASVsum >= minTaxaSum) >= minNoSubj, 
  prevalent = colSums(ASVprev >= minTaxaPrev) >= minNoSubj) %>%
  filter(abundant, prevalent)

diet <- prune_taxa(keepASV$Seq_ID, diet)
diet
diet_SMP <- data.frame(diet@sam_data) %>%
  mutate(Samp_Date = as.Date(Samp_Date))

sample_data(diet)$type <- with(
  sample_data(diet), 
  memisc::cases(
    Diet_RelDay < -15                    -> "Ancient",
    Diet_RelDay < 0 & Diet_RelDay >= -15 -> "PreDiet",
    Diet_RelDay < 5 & Diet_RelDay >= 0   -> "MidDiet",
    Diet_RelDay < 15 & Diet_RelDay >= 5  -> "PostDiet",
    Diet_RelDay >= 15                    -> "Recovery"
  )
)

## Abx Samples ----------------------------------------------------------------

abx <- subset_samples(ps_norm, grepl("Abx", Group))
abx <- subset_samples(
  abx, 
  !grepl("MidDiet", Diet_Interval) & !grepl("PreDiet", Diet_Interval) &
  !grepl("PreCC", CC_Interval))

# Filter taxa that occurs at minTaxaSum in at least minNoSubj subjects
minTaxaSum <- 10; minTaxaPrev <- 3; minNoSubj <- 3
ASVsum <- data.frame(t(as(otu_table(abx), "matrix"))) %>%
  mutate(Subject = abx@sam_data$Subject) %>%
  group_by(Subject) %>%
  summarise_all(sum)
ASVprev <- data.frame(t(as(otu_table(abx), "matrix")) > 0) %>%
  mutate(Subject = abx@sam_data$Subject) %>%
  group_by(Subject) %>%
  summarise_all(sum)

keepASV <- data.frame(
  Seq_ID = colnames(ASVsum),
  abundant = colSums(ASVsum >= minTaxaSum) >= minNoSubj, 
  prevalent = colSums(ASVprev >= minTaxaPrev) >= minNoSubj) %>%
  filter(abundant, prevalent)

abx <- prune_taxa(keepASV$Seq_ID, abx)
abx
abx_SMP <- data.frame(abx@sam_data) %>%
  mutate(Samp_Date = as.Date(Samp_Date))

sample_data(abx)$type <- with(
  sample_data(abx), 
    memisc::cases(
      Abx_RelDay < -20                   -> "Ancient",
      Abx_RelDay < 0 & Abx_RelDay >= -20 -> "PreAbx",
      Abx_RelDay < 7 & Abx_RelDay >= 0   -> "MidAbx",
      Abx_RelDay < 30 & Abx_RelDay >= 7  -> "PostAbx",
      Abx_RelDay >= 30                   -> "Recovery"
    )
)

################################################################################
## Tree Discriminant Analysis
################################################################################

## Two classes only ============================================================

## Response to Diet -----------------------------------------------------------

### Two classes only
diet_change <- subset_samples(diet, type %in%  c("PostDiet", "PreDiet")) 
diet_change <- subset_taxa(diet_change, taxa_sums(diet_change) > 0)
diet_change
cat("Starting diet two classes treeDA...")
set.seed(0)
diet_change.treedacv = treedacv(
  response = sample_data(diet_change)$type,
  predictors = asinh(t(otu_table(diet_change))),
  tree = phy_tree(diet_change),
  folds = 5, pvec = seq(50, 150, by = 10))
cat("Runtime: ", proc.time() - t, "\n")
diet_change.treedacv

save(list = c("diet_change.treedacv", "diet_change"), file = resfile)

## Response to Antibiotics -----------------------------------------------------

abx_change <- subset_samples(abx, type %in%  c("PostAbx", "PreAbx")) 
abx_change <- subset_taxa(abx_change, taxa_sums(abx_change) > 0)
abx_change

## Cross Validation
cat("Starting abx two classes treeDA...")
set.seed(0)
abx_change.treedacv = treedacv(
  response = sample_data(abx_change)$type,
  predictors = asinh(t(otu_table(abx_change))),
  tree = phy_tree(abx_change),
  folds = 5, pvec = seq(50, 150, by = 10))
cat("Runtime: ", proc.time() - t, "\n")
abx_change.treedacv

save(list = c("diet_change.treedacv", "diet_change",
              "abx_change.treedacv", "abx_change"), 
     file = resfile)

## Multiple classes ============================================================

## Response to Diet ------------------------------------------------------------

cat("Starting diet multiple treeDA...\n")
t <- proc.time()
set.seed(0)
diet.treedacv = treedacv(
  response = sample_data(diet)$type,
  predictors = asinh(t(otu_table(diet))),
  tree = phy_tree(diet),
  folds = 5, pvec = seq(50, 150, by = 10))
cat("Runtime: ", proc.time() - t, "\n")
diet.treedacv

save(list = c("diet.treedacv", "diet", "diet_change.treedacv", "diet_change",
              "abx_change.treedacv", "abx_change"), 
     file = resfile)

## Response to Antibiotics -----------------------------------------------------

cat("Starting abx multiple treeDA...")
set.seed(0)
abx.treedacv = treedacv(
  response = sample_data(abx)$type,
  predictors = asinh(t(otu_table(abx))),
  tree = phy_tree(diet),
  folds = 5, pvec = seq(50, 150, by = 10))
cat("Runtime: ", proc.time() - t, "\n")
abx.treedacv

save(list = c("diet.treedacv", "diet", "diet_change.treedacv", "diet_change",
              "abx.treedacv", "abx", "abx_change.treedacv", "abx_change"), 
     file = resfile)

