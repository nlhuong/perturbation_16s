#!/usr/bin/env Rscript
##
## File description -------------------------------------------------------------
##
## Generate a phylogenetic tree for sequences in the phyloseq object.
##
## author: nlhuong90@gmail.com
## date: 12/21/2017

library(phyloseq)
library(DECIPHER)
library(phangorn)

ALIGN <- FALSE
DIST <- FALSE
PHYLONJ <- FALSE
PHYLOGTR <- TRUE

path_to_process_data <- "/home/lanhuong/Projects/PerturbationStudy/data/processed"
data_file <- "perturb_physeq_filtered_8Dec.rds"
res_filename <- "phanghorn_fit.rda"
res_prev_run <- "phanghorn_fit_0.rda"
res_alignfile <- "decipher_align_0.rds"

ps <- readRDS(file.path(path_to_process_data, data_file))
seqs <- as.character(ps@tax_table[, "Seq"])
names(seqs) <- taxa_names(ps)
head(seqs)

if(ALIGN) {
  ptm <- proc.time()
  alignment <- AlignSeqs(DNAStringSet(seqs), anchor= NA, verbose=FALSE)
  alg.time <- proc.time() - ptm
  saveRDS(alignment, file.path(path_to_process_data, res_alignfile))
  cat("Completed alingment in: \n")
  print(alg.time)
} else {
  alignment <- readRDS(file.path(path_to_process_data, res_alignfile))
  cat("Loaded alingment. \n")
}

if(DIST) {
  ptm0 <- proc.time()
  phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
  dm <- dist.ml(phang.align)
  cat("Computed dist ml.n:\n")
  print(proc.time() - ptm0)
  ptm <- proc.time()
  save(list = c("alignment", "phang.align", "dm"), 
       file = file.path(path_to_process_data, res_filename))
}

if(PHYLONJ) {
  treeNJ <- NJ(dm) 
  fit <- pml(treeNJ, data = phang.align)
  cat("Computed NJ tree in:\n")
  print(proc.time() - ptm)
  ptm <- proc.time()
  save(list = c("alignment", "phang.align", "dm", "treeNJ", "fit"), 
       file = file.path(path_to_process_data, res_filename))
} 

if(!all(c(DIST, PHYLONJ))) {
  load(file.path(path_to_process_data, res_prev_run))
  cat("Loaded objects. \n")
}

if(PHYLOGTR){
  ptm <- proc.time()
  fitGTR <- update(fit, k=4, inv=0.2)
  fitGTR <- optim.pml(fitGTR, model="GTR", optInv = TRUE, optGamma = TRUE,
                      rearrangement = "stochastic") 
  cat("Computed GTR tree in:\n")
  print(proc.time() - ptm)
  save(list = c("alignment", "phang.align", "dm", "treeNJ", "fit", "fitGTR"), 
       file = file.path(path_to_process_data, res_filename))
}
