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
load("results/ordinate_16S.rda")


# Tree Discriminant Analysis

## Response to Diet

dds <- phyloseq_to_deseq2(ps, design = ~ 1)
dds <- estimateSizeFactors(dds, type = "poscounts")
norm_counts <- counts(dds, normalized = TRUE)
ps_norm <- ps
otu_table(ps_norm) <- otu_table(norm_counts, taxa_are_rows = TRUE)


minTaxaSum <- 5
diet <- subset_samples(ps_norm, grepl("Diet", Group))
diet <- subset_samples(
  diet, !grepl("MidAbx", Abx_Interval) & !grepl("PostAbx", Abx_Interval) &
    !grepl("PostCC", CC_Interval))

diet <- subset_taxa(diet, taxa_sums(diet) > minTaxaSum)
diet_SMP <- data.frame(diet@sam_data) %>%
  mutate(Samp_Date = as.Date(Samp_Date))


sample_data(diet)$type <- with(
  sample_data(diet), 
  memisc::cases(
     Diet_RelDay < 0                      -> "PreDiet",
     Diet_RelDay < 7 & Diet_RelDay >= 0   -> "MidDiet",
     Diet_RelDay < 40 & Diet_RelDay >= 7  -> "PostDiet",
     Diet_RelDay >= 40                    -> "Recovery"
   )
)


# Cross Validation
cat("Starting diet multiple treeDA...\n")
t <- proc.time()
set.seed(0)
diet.treedacv = treedacv(
  response = sample_data(diet)$type,
  predictors = t(otu_table(diet)),
  tree = phy_tree(diet),
  folds = 5, pvec = seq(20, 200, by = 20))
cat("Runtime: ", proc.time() - t, "\n")
diet.treedacv

save(list = c("diet.treedacv", "diet"), file = resfile)

### Two classes only
  
diet_change <- subset_samples(diet, type %in%  c("PostDiet", "PreDiet")) 
diet_change <- subset_taxa(diet_change, taxa_sums(diet_change) > 0)
diet_change
#  Cross Validation
cat("Starting diet two classes treeDA...")
set.seed(0)
diet_change.treedacv = treedacv(
  response = sample_data(diet_change)$type,
  predictors = t(otu_table(diet_change)),
  tree = phy_tree(diet_change),
  folds = 5, pvec = seq(20, 200, by = 20))
cat("Runtime: ", proc.time() - t, "\n")
diet_change.treedacv

save(list = c("diet.treedacv", "diet", "diet_change.treedacv", "diet_change"), file = resfile)

## Response to Antibiotics
    
minTaxaSum <- 5
abx <- subset_samples(ps_norm, grepl("Abx", Group))
abx <- subset_samples(abx, 
                      !grepl("MidDiet", Diet_Interval) & !grepl("PreDiet", Diet_Interval) &
                        !grepl("PreCC", CC_Interval))

abx <- subset_taxa(abx, taxa_sums(abx) > minTaxaSum)
abx_SMP <- data.frame(abx@sam_data) %>%
  mutate(Samp_Date = as.Date(Samp_Date))
abx


### Multiple Classes
    
sample_data(abx)$type <- with(
  sample_data(abx), 
    memisc::cases(
      Abx_RelDay < 0                    -> "PreAbx",
      Abx_RelDay < 7 & Abx_RelDay >=0   -> "MidAbx",
      Abx_RelDay < 40 & Abx_RelDay >= 7 -> "PostAbx",
      Abx_RelDay >= 40                  -> "Recovery"
    )
)


##  Cross Validation
cat("Starting abx multiple treeDA...")
set.seed(0)
abx.treedacv = treedacv(
  response = sample_data(abx)$type,
  predictors = t(otu_table(abx)),
  tree = phy_tree(diet),
  folds = 5, pvec = seq(20, 200, by = 20))
cat("Runtime: ", proc.time() - t, "\n")
abx.treedacv

save(list = c("diet.treedacv", "diet", "diet_change.treedacv", "diet_change",
              "abx.treedacv", "abx"), 
file = resfile)

      
      
### Two classes only
      
abx_change <- subset_samples(abx, type %in%  c("PostAbx", "PreAbx")) 
abx_change <- subset_taxa(abx_change, taxa_sums(abx_change) > 0)
abx_change

## Cross Validation
cat("Starting abx two classes treeDA...")
set.seed(0)
abx_change.treedacv = treedacv(
  response = sample_data(abx_change)$type,
  predictors = t(otu_table(abx_change)),
  tree = phy_tree(abx_change),
  folds = 10, pvec = seq(20, 200, by = 20))
cat("Runtime: ", proc.time() - t, "\n")
abx_change.treedacv

save(list = c("diet.treedacv", "diet", "diet_change.treedacv", "diet_change",
                      "abx.treedacv", "abx", "abx_change.treedacv", "abx_change"), 
             file = resfile)

        
