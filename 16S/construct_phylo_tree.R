#!/usr/bin/env Rscript
##
## File description -------------------------------------------------------------
##
## Generate a phylogenetic tree for sequences in the phyloseq object.
##
## author: nlhuong90@gmail.com
## date: 12/21/2017; updated 07/11/2018

library(phyloseq)
library(DECIPHER)
library(phangorn)

ALIGN <- FALSE
DIST <- FALSE
PHYLONJ <- FALSE
PHYLOGTR <- TRUE

path_to_process_data <- "/scratch/users/lanhuong/Projects/perturbation_16s/data/16S/phyloseq/"
datafile <- file.path(path_to_process_data, "perturb_physeq_fltr2438_15Jul19.rds") 
resfile <- file.path(path_to_process_data, "phylo_fit_fltr2438_15Jul19.rda")


if(file.exists(resfile)) {
  load(resfile)
}

ps <- readRDS(datafile)
ps
seqs <- as.character(ps@tax_table[, "Seq"])
names(seqs) <- taxa_names(ps)
head(seqs)

if(ALIGN & !("alignment" %in% ls())) {
  ptm <- proc.time()
  alignment <- AlignSeqs(DNAStringSet(seqs), anchor= NA, verbose=FALSE)
  alignment.time <- proc.time() - ptm
  save(list = c("alignment", "alignment.time"), file = resfile)
  cat("Completed alingment in: \n")
  print(alignment.time)
} 



if(DIST & !("dm" %in% ls())) {
  ptm0 <- proc.time()
  phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
  dm <- dist.ml(phang.align)
  cat("Computed dist ml.n:\n")
  dist.time <- proc.time() - ptm0
  print(dist.time)
  save(list = c("alignment", "alignment.time",
                "phang.align", "dm", "dist.time"), 
       file = resfile)
}

if(PHYLONJ & !("treeNJ" %in% ls())) {
  ptm0 <- proc.time()
  treeNJ <- NJ(dm) 
  fitNJ <- pml(treeNJ, data = phang.align)
  rootedNJtree <- phangorn::midpoint(fitNJ$tree)
  cat("Computed NJ tree in:\n")
  treenj.time <- proc.time() - ptm0
  print(treenj.time)
  save(list = c("alignment", "alignment.time",
                "phang.align", "dm", "dist.time", 
                "treeNJ", "fitNJ", "rootedNJtree", "treenj.time"), 
       file = resfile)
} 


if(PHYLOGTR  & !("fitGTR" %in% ls())){
  ptm0 <- proc.time()
  fitGTR <- update(fitNJ, k=4, inv=0.2)
  fitGTR <- optim.pml(fitGTR, model="GTR", optInv = TRUE, optGamma = TRUE,
                      rearrangement = "stochastic") 
  rootedGTRtree <- phangorn::midpoint(fitGTR$tree)
  cat("Computed GTR tree in:\n")
  treegtr.time <- proc.time() - ptm0
  print(treegtr.time)
  save(list = c("alignment", "alignment.time",
                "phang.align", "dm", "dist.time", 
                "treeNJ", "fitNJ", "rootedNJtree", "treenj.time", 
                "fitGTR", "rootedGTRtree", "treegtr.time"), 
       file = resfile)
}




