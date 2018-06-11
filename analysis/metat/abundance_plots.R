#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Plotting counts. 
##
## author: nlhuong90@gmail.com
## data 3/15/201

###############################################################################
## Libraries and setup
###############################################################################
library("tidyverse")
library("readxl")
library("data.table")

datadir <- "/scratch/PI/sph/resilience/metatranscriptomics/processed/"
setwd(datadir)

theme_set(theme_bw())
theme_updates <- theme(text = element_text(size =20))

###############################################################################
## Load data
###############################################################################

metat <- fread("final_results/abund_refseq_len_ratio.csv")
gene_info <- fread("final_results/gene_len_refseq_raw.csv")

dim(metat)
#[1] 1278950     383
dim(gene_info)
#[1] 1278950       4


gene_info2 <- metat[, .(GeneID, Organism, Function)]
identical(gene_info2[order(GeneID)], 
          gene_info[order(GeneID), !c("Length"), with = FALSE])
rm(gene_info2)
metat <- metat[, !c("Organism", "Function"), with = FALSE]

gene_stats <- metat[, 
                    .(GeneID, 
                      sum = rowSums(.SD, na.rm = TRUE),
                      prev = rowSums(.SD > 0, na.rm = TRUE)
                      ), 
                    .SDcols = setdiff(colnames(metat), "GeneID")] 

summary(gene_stats)
# GeneID               sum                prev        
# Length:1278950     Min.   :    0.00   Min.   :  1.000  
# Class :character   1st Qu.:    0.06   1st Qu.:  1.000  
# Mode  :character   Median :    0.13   Median :  2.000  
                    # Mean   :    1.13   Mean   :  6.963
                    # 3rd Qu.:    0.39   3rd Qu.:  5.000
                    # Max.   :37487.68   Max.   :379.000
identical(gene_stats$GeneID, metat$GeneID)
metat <- metat[gene_stats$prev >1, ]
dim(metat)
#[1] 746555    381

metat <- metat[GeneID %in% gene_stats[prev > 5, GeneID], ]
dim(metat)
#[1] 318372    381

###############################################################################
## Load data
###############################################################################























