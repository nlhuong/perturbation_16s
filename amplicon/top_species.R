#!/usr/bin/env Rscript
##
## File description -------------------------------------------------------------
##
## List top species comprising a fix fraction of the total number of reads
## within a sample or collections of samples (e.g. from same subject).
##
## author: nlhuong90@gmail.com
## date: 01/15/2018

###############################################################################
## Setup and data
###############################################################################
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(viridis)

options(stringsAsFactors = FALSE)
theme_set(theme_bw())
theme_update(text = element_text(size = 15))

path2data <- "perturbation_16s/data/processed/"
path2figs <- "perturbation_16s/figs"
path2out <- "perturbation_16s/output"

ps0 <- ps <- readRDS(file.path(path2data, "perturb_physeq_filtered_27Dec.rds"))

## Function that returns for each sample (column), the top sequences (rows) 
## covering a specified fixed fraction of all reads in the sample. 
top_seqs <- function(seqtab, frac){
  apply(seqtab, 2, function(col) {
    ord_idx <- order(col, decreasing = TRUE)
    csum <- cumsum(col[ord_idx])
    ntaxa <- sum(csum < frac * sum(col)) 
    return(rownames(seqtab)[ord_idx[1:ntaxa]])
  })
}

extract_taxa <- function(seqs, taxtab, taxrank = "Species") {
  colnames(taxtab) <- tolower(colnames(taxtab))
  col_idx <- grep(tolower(taxrank), colnames(taxtab))
  taxa <- taxtab[seqs, col_idx]
  # Find non-NA taxrank related names corresponding to seqs in taxtab
  taxa <- apply(taxa, 1, function(x) {x[which(!is.na(x))[1]]})
  return(taxa)
}

get_stats <- function(toptaxa, sample_data) {
  nseqs <- sapply(toptaxa[["seqs"]], length)
  nspecies <- sapply(toptaxa[["species"]], 
                     function(x) length(unique(x[!is.na(x)])))
  ngenus <- sapply(toptaxa[["genus"]], 
                   function(x) length(unique(x[!is.na(x)])))
  data.frame(nseqs, nspecies, ngenus) %>%
    rownames_to_column("Subject") %>%
    arrange(nseqs) %>%
    left_join(sample_data)
}

###############################################################################
## # (Optional) agglomerate samples by criteria from sample data
###############################################################################

# Here we agglomerate samples from same subjects
groups <- "Subject"
sample_data <- data.frame(ps0@sam_data)

group_data <- sample_data %>%
  mutate(depth = sample_sums(ps0)) %>%
  select(Subject, Group, depth) %>%
  group_by(Subject, Group) %>%
  summarise(
    mean_sample_sum = mean(depth),
    total_read_sum  =sum(depth)
  ) %>%
  arrange(Group) %>%
  as.data.frame() %>%
  column_to_rownames("Subject") 
  
group_seqtab <- sapply(unique(sample_data[[groups]]), function(g) 
  rowSums(ps@otu_table[, sample_data[[groups]] == g]))

group_seqtab <- group_seqtab[, rownames(group_data)]
ps <- phyloseq(otu_table(group_seqtab, taxa_are_rows = TRUE),
               sample_data(group_data), tax_table(ps))

###############################################################################
## Find top taxa covering fixed fraction of reads
###############################################################################

# coverage fraction
frac <- 0.95

seqtab <- as(ps@otu_table, "matrix")
taxtab <- data.frame(ps@tax_table)
  
toptaxa <- list(seqs = top_seqs(seqtab, frac))
toptaxa[["species"]] <- lapply(toptaxa$seqs, function(seqs) 
    extract_taxa(seqs, taxtab, "species"))
toptaxa[["genus"]] <- lapply(toptaxa$seqs, function(seqs) 
    extract_taxa(seqs, taxtab, "genus"))

toptaxa_df <- lapply(toptaxa, function(lst) {
  max_len <- max(sapply(lst, length))
  lst <- lapply(lst, function(x) {
    x[is.na(x)] <- "Unknown"
    c(x, rep(NA, max_len- length(x)))
  })
  return(data.frame(lst))
})

write.csv(toptaxa_df$seqs, file.path(path2out, paste0("seqs_", 100*frac,".csv")))
write.csv(toptaxa_df$species, file.path(path2out, paste0("species_", 100*frac,".csv")))
write.csv(toptaxa_df$genus, file.path(path2out, paste0("genus_", 100*frac,".csv")))

###############################################################################
## Visualize the data and save results
###############################################################################

toptaxa_stat <- get_stats(toptaxa, group_data %>% rownames_to_column("Subject"))
plt <- ggplot(toptaxa_stat %>% filter(!Subject == "NA_QC")) +
  geom_text(
    size = 5,
    aes(x = log10(mean_sample_sum), 
        y = nseqs,
        color = log10(total_read_sum), 
        label = Subject)) +
  scale_color_viridis(end = 0.97) + 
  ylab(paste0("No. of seqs covering ", 100*frac, "% reads"))

filename <- paste0("toptaxa_", 100*frac, "coverage.png")
ggsave(plt, filename = file.path(path2figs, filename))

