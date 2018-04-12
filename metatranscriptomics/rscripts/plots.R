#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Plotting counts. 
##
## author: nlhuong90@gmail.com
## date: 03/15/2018

###############################################################################
## Libraries and setup
###############################################################################
library("tidyverse")
library("readxl")
library("phyloseq")

theme_set(theme_bw())
theme_updates <- theme(text = element_text(size =20))

rpkm <- function(counts, lengths) {
  rate <- counts / lengths 
  rate / sum(counts) * 1e6
}

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

plt_trend <- function(df, x, y){
  dfrepel <- df
  if(nrow(df) > 10^4) {
    dfrepel <- df[sample(1:nrow(df), floor(nrow(df)/2)), ]
  }
  ggplot(
    df, aes_string(x, y)
  ) +
    stat_density_2d(
      aes(fill = ..density..), 
      geom = "raster", contour = FALSE) +
    geom_point(
      data = dfrepel,
      alpha = 0.1, size = 0.7) +
    scale_fill_distiller(palette= "Spectral", direction=-1) +
    scale_x_continuous(
      expand = c(0, 0), trans = "log10") +
    scale_y_continuous(
      limits = c(NA, quantile(df[[y]], 0.95)),
      expand = c(0, 0)) 
}

plot_abundance <- function(X, figfile, title, N = 10^5) {
  geneStats <- data.frame(
    GeneID = rownames(X),
    GeneSum = rowSums(X),
    GeneMean = rowMeans(X),
    GenePrev = rowSums(X > 0)/ncol(X),
    GeneSD = apply(X, 1, sd)
  )
  dfplot <- geneStats %>% 
    filter(GeneMean < quantile(geneStats$GeneMean, 0.99))
  if(nrow(X) > N) {
    dfplot <- dfplot[sample(1:nrow(dfplot), N), ]
  }
  # plots
  # Sample Depth plot
  plt1 <- qplot(colSums(X)) +
    ggtitle(paste0("[", title,  "] sample depths"))
  # Mean vs revalence plot
  plt2 <- plt_trend(dfplot, "GeneMean", "GenePrev") + 
    ggtitle(paste0("[", title,  "] mean vs prevalence"))
  # Mean vs variance plot
  plt3 <- plt_trend(dfplot, "GeneMean", "GeneSD") + 
    ggtitle(paste0("[", title,  "] mean vs sd"))
  pdf(figfile, width = 8, height = 5)
  print(plt1)
  print(plt2)
  print(plt3)
  dev.off()
}

################################################################################

meas <- read_xlsx("data/Mapping_Files_22Jan2018.xlsx", "Meas", skip = 1) %>%
  rename(Samp_ID = SampID)
interv_levs <- c("NoInterv", "PreDiet", "MidDiet", "PostDiet", "PreCC", "MidCC",
                 "PostCC", "PreAbx", "MidAbx", "PostAbx")
samp <- read_xlsx("data/Mapping_Files_22Jan2018.xlsx", "Samp", skip = 1) %>%
  mutate(
    Diet_Interval = ifelse(Diet_Interval == "NA", "NoInterv", Diet_Interval),
    CC_Interval = ifelse(CC_Interval == "NA", "NoInterv", CC_Interval),
    Abx_Interval = ifelse(Abx_Interval == "NA", "NoInterv", Abx_Interval),
    Diet_Interval = factor(Diet_Interval, interv_levs),
    CC_Interval = factor(CC_Interval, interv_levs),
    Abx_Interval = factor(Abx_Interval, interv_levs)
  ) %>%
  filter(Samp_Type != "ExtrCont")

meas <- meas %>% select(Meas_ID, Samp_ID) %>%
  left_join(samp)
rm(samp)

# loading counts
# MetatT Data
metat <- read_csv("data/metatranscriptomic/counts/Annotated_Relman_RNAseq_16_refseq_genes.csv") 
# 1774974     100
prevalence <- rowSums(metat %>% select(starts_with("M")) > 0)
metat <- metat[prevalence > 5, ]
# 606694    100
geneSums <- metat %>% select(starts_with("M")) %>% rowSums()
metat <- metat[geneSums > 10, ]
# 491124    100
smp_names <- metat %>%select(starts_with("M")) %>% colnames() 
meas_names <- gsub("\\_(.*)", "", smp_names)
meas <- meas %>% filter(Meas_ID %in% meas_names)
table(meas$Subject)
# DBU EBF 
# 10  84 

curr_subj <- "EBF"
figfile <- paste0("figs/metatrans/", curr_subj)

meas <- meas %>% filter(Subject == curr_subj)

metatInfo <- metat %>%
  select(GeneID, Length, Organism, Function) %>%
  separate(Organism, c("Genus", "species"), sep = " ", remove = FALSE,
           extra = "merge", fill = "right") %>%
  rename(species_id = Organism)


RAW <- metat %>%select(starts_with("M")) %>% as.matrix()
rownames(RAW) <- metatInfo$GeneID
colnames(RAW) <- gsub("\\_(.*)", "", smp_names)
RAW <- RAW[, meas$Meas_ID]
RAW <- RAW[rowSums(RAW) > 1, ]
metatInfo <- metatInfo %>% filter(GeneID %in% rownames(RAW))

rm(prevalence, geneSums, metat)


################################################################################

plot_abundance (RAW, paste0(figfile, "_raw_abund.pdf"), 
                paste0("Subj. ", curr_subj, ": RAW Counts"))
  
################################################################################


RPKM <- apply(RAW, 2, function(x) rpkm(x, metatInfo$Length))
plot_abundance (RPKM, paste0(figfile, "_rpkm_abund.pdf"),
                paste0("Subj. ", curr_subj, ": RPKM Counts"))
rm(RPKM)

################################################################################

TPM <- apply(RAW, 2, function(x) tpm(x, metatInfo$Length))
plot_abundance (TPM, paste0(figfile, "_tpm_abund.pdf"),
                paste0("Subj. ", curr_subj, ": TPM Counts"))

################################################################################

ps <- readRDS("data/16s/perturb_physeq_22Jan18.rds") 
taxtab <- data.frame(ps@tax_table) 
taxtab <- taxtab %>% select(-Seq)
rm(ps)

org_tpm <- aggregate(TPM, by = list(metatInfo$species_id), FUN = sum)
org_tpm <- org_tpm %>%  column_to_rownames("Group.1")
org_tpm_sums <- rowSums(org_tpm)

plot_abundance(org_tpm, paste0(figfile, "_tpm_species_abund.pdf"), 
                paste0("Subj. ", curr_subj, ": species TPM species"))

idx <- order(-org_tpm_sums)
org_tpm <- org_tpm[idx, ]
interval <- pretty(1:length(org_tpm_sums), 6)
keep_idx <- c(sample(1:interval[2], 10), 
              sample(interval[3]:interval[4], 10),
              sample(interval[5]:nrow(org_tpm), 10))
X <- org_tpm[keep_idx, ]

melt_org_tmp <- X %>%
  rownames_to_column("species_id") %>%
  gather(Meas_ID, TPM, -species_id) %>%
  separate(species_id, c("Genus", "species"), 
           sep = " ", remove = FALSE,
           extra = "merge", fill = "right") %>%
  left_join(meas %>% select("Meas_ID", "Samp_ID", "Subject",
                            "Abx_Interval", "Diet_Interval")) %>%
  mutate(species_id = factor(species_id, levels = rownames(X)))

pdf("figs/metatrans/EBF_orgs_tpm.pdf", width = 8, height = 6)
ggplot(melt_org_tmp) +
  geom_histogram(
    aes(x = TPM),
    bins = 30) +
  facet_wrap(~ species_id, scale = "free", ncol = 5)
dev.off()


genus_tpm <- aggregate(TPM, by = list(metatInfo$Genus), FUN = sum)
genus_tpm <- genus_tpm %>%  column_to_rownames("Group.1")
genus_tpm_sums <- rowSums(genus_tpm)

idx <- order(-genus_tpm_sums, na.last = NA)
genus_tpm <- genus_tpm[idx, ]
interval <- seq(1, length(genus_tpm_sums), length.out = 5)
keep_idx <- c(1:10, interval[3]:(interval[3]+9), interval[4]:(interval[4] + 9))
X <- genus_tpm[keep_idx, ]

genustab <- taxtab[, 1:6]
genustab <- genustab[!duplicated(genustab), ]

melt_genus_tmp <- X %>%
  rownames_to_column("Genus") %>%
  gather(Meas_ID, TPM, -Genus) %>%
  left_join(genustab) %>%
  mutate(Genus = factor(Genus, levels = rownames(X)))

pdf("figs/metatrans/EBF_genus_tpm.pdf", width = 8, height = 6)
ggplot(melt_genus_tmp) +
  geom_histogram(
    aes(x = TPM, fill = Family),
    bins = 30) +
  facet_wrap(~ Genus, scale = "free", ncol = 5)
dev.off()



################################################################################
function_tpm <- aggregate(TPM, by = list(metatInfo$Function), FUN = sum)
function_tpm <- function_tpm %>% column_to_rownames("Group.1")

fun_sums <-  rowSums(function_tpm)
idx <- order(-fun_sums)
function_tpm <- function_tpm[idx, ]
interval <- pretty(1:length(fun_sums), 6)
keep_idx <- c(sample(1:interval[2], 10), 
              sample(interval[3]:interval[4], 10),
              sample(interval[5]:nrow(function_tpm), 10))
X <- function_tpm[keep_idx, ]

melt_fun_tmp <- X %>%
  rownames_to_column("function_id") %>%
  gather(Meas_ID, TPM, -function_id) %>%
  left_join(meas %>% select("Meas_ID", "Samp_ID", "Subject",
                            "Abx_Interval", "Diet_Interval")) %>%
  mutate(function_id = factor(function_id, levels = rownames(X)))

pdf("figs/metatrans/EBF_fun_tpm.pdf", width = 8, height = 6)
ggplot(melt_fun_tmp) +
  geom_histogram(
    aes(x = TPM),
    bins = 30) +
  facet_wrap(~ function_id, scale = "free", ncol = 5)+
  theme(text = element_text(15))
function_iddev.off()



library("tidyverse")
metatInfoSub <- metatInfo %>%
  filter(GeneID %in% rownames(TPM))


# RAW counts
# summary(geneStats)
# GeneID          GeneSum            GeneMean           GenePrev         GeneSD        
# WP_000001645.1:     1   Min.   :     2.0   Min.   :   0.024   Min.   : 1.00   Min.   :   0.153  
# WP_000003861.1:     1   1st Qu.:    15.0   1st Qu.:   0.179   1st Qu.:11.00   1st Qu.:   0.526  
# WP_000020194.1:     1   Median :    29.0   Median :   0.345   Median :18.00   Median :   0.824  
# WP_000029030.1:     1   Mean   :   147.4   Mean   :   1.755   Mean   :24.51   Mean   :   2.778  
# WP_000029031.1:     1   3rd Qu.:    76.0   3rd Qu.:   0.905   3rd Qu.:33.00   3rd Qu.:   1.682  
# WP_000029276.1:     1   Max.   :660890.0   Max.   :7867.738   Max.   :84.00   Max.   :4337.447  
# (Other)       :468429v







