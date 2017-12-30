
## File description -------------------------------------------------------------
##
## Generate reads quality plots (QC plots) for each forward and reverse files 
## in each pool in the specified 'path.split.libraries'
##
## author: nlhuong90@gmail.com
## date: 12/21/2017
## 
## Run this command on the icme-share server
# { time Rscript --no-save --no-restore --verbose dada2_qc_plots.R > 
# qc_plots.rout 2> qc_plots.err ; } 2> qc_plots.time

library(dada2)
library(doParallel)
library(foreach)
library(ggplot2)



path.split.libraries <- 
  "/home/lanhuong/Projects/PerturbationStudy/data/HMD_16S/split_libraries"

path.res <- "/home/lanhuong/Projects/PerturbationStudy/output/16S_qc_plots"
dir.create(path.res, showWarnings = FALSE)
pool.list <- c("16S1_idemp", "16S2_idemp", "16S3_idemp", "16S4_idemp")

# Register cluster
registerDoParallel(length(pool.list))
foreach(pool = pool.list) %dopar% {
  # Sort ensures forward/reverse reads are in same order
  pDir <- file.path(path.split.libraries, pool)
  pool.name <- gsub("_idemp", "", pool)
  fwdFQ <- sort(list.files(file.path(pDir, paste0(pool.name, "_F"))))
  revFQ <- sort(list.files(file.path(pDir, paste0(pool.name, "_R"))))
  # Specify the full path to the fwdFQ and revFQ
  fwdFQ <- file.path(pDir, paste0(pool.name, "_F"), fwdFQ)
  revFQ <- file.path(pDir, paste0(pool.name, "_R"), revFQ)
  
  pdf(file = file.path(path.res, paste0(pool.name, "_fwd_qc.pdf")))
  for (i in 1:length(fwdFQ)) {
    print(plotQualityProfile(fwdFQ[i]) + ylim(0, 42))
  }
  dev.off()
  
  pdf(file = file.path(path.res, paste0(pool.name, "_rev_qc.pdf")))
  for (i in 1:length(revFQ)) {
    print(plotQualityProfile(revFQ[i]) + ylim(0, 42))
  }
  dev.off()
}



