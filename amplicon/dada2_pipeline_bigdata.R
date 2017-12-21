
## File description -------------------------------------------------------------
##
## Use DADA2 pipeline to generate a sequence abundance table and taxonomy
## assignment from reads fastq files.
##
## author: nlhuong90@gmail.com
## date: 12/21/2017
##
## Run this command on the icme-share server
# { time Rscript --no-save --no-restore --verbose dada2_pipeline_bigdata.R > 
## dada2.rout 2> dada2.err ; } 2> dada2.time

pool.list <- c("16S1_idemp", "16S2_idemp", "16S3_idemp","16S4_idemp")

FILTER <- FALSE
LEARN_ERRORS <- FALSE
INFER_SEQUENCE_VARIANTS <- FALSE
MAKE_SEQTAB <- FALSE
COMBINE <- FALSE
CHIMERA_CHECK <- FALSE
ASSIGN_TAXONOMY <- TRUE

nCores <- 60 #for chimera checking

## Load Package
.packages <- c("dada2", "doParallel", "foreach", "ggplot2", "dplyr", "phyloseq")
sapply(.packages, require, character.only = TRUE)

## SETTINGS ---------------------------
# DADA2 Paramers
## Filter and trim parameters
seed <- 74658
maxEE <- c(2, 2)
trimLeft <- c(5, 5) 
learnErrorsNreads <- 2e+06

# Change the two below:
ncharTabTrim <- seq(350, 370)
truncLen.lst <- list("16S1_idemp" = c(240, 200), 
                     "16S2_idemp" = c(240, 200), 
                     "16S3_idemp" = c(220, 180),
                     "16S4_idemp" = c(220, 180))

## Main directory to 16S data
main.dir <- path <- "/home/lanhuong/Projects/PerturbationStudy/data/HMD_16S"

## Directory to dada2 results
path.to.res <- "/home/lanhuong/Projects/PerturbationStudy/output/16S_dada2"
dir.create(path.to.res, showWarnings = FALSE)

## Directory to reference fasta files
path.to.reference <- file.path(main.dir, "reference")

## Directory for trimmed files
path.to.split.lib <- file.path(main.dir, "split_libraries")

## Directory for trimmed files
path.to.filt <- file.path(main.dir, "filtered_libraries")
dir.create(path.to.filt, showWarnings = FALSE)


## Save all parameters to file in results directory
paramaters <- c("DADA2 parameters different from the default:", 
                paste0("[filterAndTrim()] maxEE: c(", maxEE[1], ",", maxEE[2], ")"),
                paste0("[filterAndTrim()] truncLengths: ", 
                       paste0(names(truncLen.lst)," = ", 
                              truncLen.lst, collapse = ", ")),
                paste0("[filterAndTrim()] trimLeft: ", trimLeft),
                paste0("[learnErrors()] nreads = ", learnErrorsNreads), 
                "[learnErrors()] randomize = TRUE",
                paste0("[makeSequenceTab()] trimming seq by length: ", min(ncharTabTrim), "-",
                       max(ncharTabTrim)))
writeLines(paramaters, con = file.path(path.to.res, "dada2_params_used.txt"))

## START DADA2 PIPELINE---------------------------"
if(FILTER) {
  registerDoParallel(length(pool.list))
  foreach(i = seq_along(pool.list)) %dopar% {
    ## Setting up file paths ---------------------------
    pool <- pool.list[i]
    pool.name <- gsub("_idemp", "", pool)
    cat("Initializing Pool: ", pool, " ---------------------------", "\n")
    path.to.filt.pool <- file.path(path.to.filt, pool)
    dir.create(path.to.filt.pool, showWarnings = FALSE)
    # Sort ensures forward/reverse reads samples are in same order
    path.to.pool.data <- file.path(path.to.split.lib, pool)
    fwdFQ <- sort(list.files(file.path(path.to.pool.data, paste0(pool.name, "_F"))))
    revFQ <- sort(list.files(file.path(path.to.pool.data, paste0(pool.name, "_R"))))
    # Extract sample names
    sample.names <- gsub("\\.(.*)", "", fwdFQ)
    # Specify the full path to the fwdFQ and revFQ
    # Specify the full path to the fwdFQ and revFQ
    fwdFQ <- file.path(path.to.pool.data, paste0(pool.name, "_F"), fwdFQ)
    revFQ <- file.path(path.to.pool.data, paste0(pool.name, "_R"), revFQ)
    # Names for the new filtered files
    filtFs <- file.path(path.to.filt.pool, paste0(sample.names, "_F_filt.fq.gz"))
    filtRs <- file.path(path.to.filt.pool, paste0(sample.names, "_R_filt.fq.gz"))
    ## Filtering and trimming ---------------------------
    cat("[", pool, "] filtering and trimming ---------------------------", "\n")
    # CONSIDER increasing maxEE if too few reads are found !!!
    out <- filterAndTrim(fwdFQ, filtFs, revFQ, filtRs, 
                         truncLen = truncLen.lst[[pool]],
                         rm.phix = TRUE,
                         trimLeft = trimLeft,    # trim left end since q usually lower
                         maxEE = maxEE,          # THIS IS DIFFERENT THAN THE DEFAULT
                         matchIDs = TRUE,        # Match ids between fwd and rev
                         multithread = TRUE)   # Run in parallel
    filtered <- data.frame(row.names = paste0(pool, "_", sample.names),
                           sample.id = paste0(pool, "_", sample.names),
                           file.name = rownames(out), out)
    # Save statistics
    write.csv(filtered, 
              file.path(path.to.filt.pool, paste0(pool, "_track_filter.csv")))
  }
}

if(INFER_SEQUENCE_VARIANTS) {
  registerDoParallel(length(pool.list))
  foreach(i = seq_along(pool.list)) %dopar% {
    pool <- pool.list[i]
    path.to.filt.pool <- file.path(path.to.filt, pool)
    # Use only the filtered samples which were created for both
    # fwd and rev, as some of the samples might be discarded after 
    # filtering & trimming step due to low quality, and low read number
    filtFs <- list.files(path.to.filt.pool, pattern = "_F_filt.fq.gz")
    filtRs <- list.files(path.to.filt.pool, pattern = "_R_filt.fq.gz")
    filtFs.smp.name <- gsub("_F_filt.fq.gz", "", filtFs)
    filtRs.smp.name <- gsub("_R_filt.fq.gz", "", filtRs)
    sample.names <- filtFs.smp.name
    if(!identical(filtFs.smp.name, filtRs.smp.name)) 
      stop("Forward and reverse files do not match.")
    filtFs <- file.path(path.to.filt.pool, paste0(sample.names, "_F_filt.fq.gz"))
    filtRs <- file.path(path.to.filt.pool, paste0(sample.names, "_R_filt.fq.gz"))
    names(filtFs) <- names(filtRs) <- sample.names

    ## Learning error rates ---------------------------
    if(LEARN_ERRORS) {
      cat("[", pool, "] learning error rates ---------------------------", "\n")
      set.seed(seed)
      errF <- learnErrors(filtFs, nreads = learnErrorsNreads, 
                          randomize = TRUE, multithread = TRUE)
      errR <- learnErrors(filtRs, nreads = learnErrorsNreads, 
                          randomize = TRUE, multithread = TRUE)
      save(list = c("errF", "errR"), 
           file = file.path(path.to.filt.pool, paste0(pool, "_err.rda")))
      # We also save plots
      pdf(file = file.path(path.to.filt.pool, paste0(pool, "_err_rates.pdf")))
      print(plotErrors(errF, nominalQ = TRUE) + 
              ggtitle(paste0(pool, ": forward Reads")))
      print(plotErrors(errR, nominalQ = TRUE) + 
              ggtitle(paste0(pool, ": reverse Reads")))
      dev.off()
    } else {
      cat("[", pool, "] loading error rates ---------------------------", "\n")
      load(file.path(path.to.filt.pool, paste0(pool, "_err.rda")))
    }
    
    ## Dereplication & Merging ---------------------------"
    cat("[", pool, "] dereplication and merging ---------------------------", "\n")
    # Sample inference and merger of paired-end reads
    mergers <- list()
    for(i in seq_along(sample.names)) {
      sam <- sample.names[i]
      cat("[", pool, "] Processing:", sam, "\n")
      derepF <- derepFastq(filtFs[[sam]])
      ddF <- dada(derepF, err=errF, multithread=TRUE)
      derepR <- derepFastq(filtRs[[sam]])
      ddR <- dada(derepR, err=errR, multithread=TRUE)
      merger <- mergePairs(ddF, derepF, ddR, derepR)
      rm(ddF, derepF, ddR, derepR)
      mergers[[paste0(pool, "_", sam)]] <- merger
    }
    saveRDS(mergers, file.path(path.to.filt.pool, paste0(pool, "_mergers.rds")))
    
    ## Track Reads ---------------------------"
    cat("[", pool, "] tracking files ---------------------------", "\n")
    # Save the reads remaining after each step 
    getN <- function(x) sum(getUniques(x))
    filtered.file <- file.path(path.to.filt.pool, paste0(pool, "_track_filter.csv"))
    filtered <- read.csv(file = filtered.file, row.names = 1)
    filtered <- filtered[as.character(names(mergers)), ]
    track <- cbind(filtered, sapply(mergers, getN))
    colnames(track) <- c("sample.id", "file.name", "input", "filtered", "merged")
    rownames(track) <- track$sample.id
    write.csv(track, file.path(path.to.filt.pool, paste0(pool, "_track.csv")))
  }
}

if (MAKE_SEQTAB) {
  cat("## Make seqtabs ---------------------------", "\n")
  registerDoParallel(length(pool.list))
  seqtab.lst <- foreach(i = seq_along(pool.list)) %dopar% {
    pool <- pool.list[[i]]
    path.to.filt.pool <- file.path(path.to.filt, pool)
    merged.reads <- file.path(path.to.filt.pool, paste0(pool, "_mergers.rds"))
    mergers <- readRDS(merged.reads)
    cat("[", pool, "] Make table ---------------------------", "\n")
    seqtab <- makeSequenceTable(mergers)
    cat("[", pool, "] Count table dimensions:", dim(seqtab), "\n")
    saveRDS(seqtab, file.path(path.to.res, paste0(pool, "_seqtab.rds"))) 
    return(seqtab)
  }
}

if(COMBINE) {
  if (length(pool.list) <= 1) 
    stop("Only one pool, nothing to combine")
  cat("## Combine seqtabs ---------------------------",  "\n")
  track.list <- list()
  seqtaball <- readRDS(file.path(path.to.res, paste0(pool.list[[1]], "_seqtab.rds")))
  for (i in 1:length(pool.list)) {
    pool <- pool.list[[i]]
    path.to.filt.pool <- file.path(path.to.filt, pool)
    track.file <- file.path(path.to.filt.pool, paste0(pool, "_track.csv"))
    track.list[[pool]] <- read.csv(track.file, row.names = 1)
    if(i != 1) {
      seqtab.file <- file.path(path.to.res, paste0(pool, "_seqtab.rds"))
      seqtab <- readRDS(seqtab.file)
      seqtaball <- mergeSequenceTables(seqtaball, seqtab)
    } 
  }
  
  cat("## Save reads tracking data ---------------------------", "\n")
  track.df <- plyr::ldply(track.list, .id = "pool")
  row.names(track.df) <- track.df$sample.id
  track.file <- file.path(path.to.res, "track_all_pools.csv")
  write.csv(track.df, file = track.file)
  
  cat("## Save sequence table ---------------------------", "\n")
  cat("Count table dimensions:", dim(seqtaball), "\n")
  saveRDS(seqtaball, file.path(path.to.res, "seqtab_all_init.rds")) 
}


if(CHIMERA_CHECK) {
  seqtaball <- readRDS(file.path(path.to.res, "seqtab_all_init.rds"))
  # Inspect distribution of sequence lengths
  cat("Initial count table dimensions:", dim(seqtaball), "\n")
  seqtaball <- seqtaball[, colSums(seqtaball) > 0]
  cat("Non-zero count table dimensions:", dim(seqtaball), "\n")
  cat("Distribution of sequence lengths: \n")
  print(table(nchar(getSequences(seqtaball))))
  seqtab.rightlen <- seqtaball[, nchar(colnames(seqtaball)) %in% ncharTabTrim]
  cat("Count table rightlen dim:", dim(seqtab.rightlen), "\n")
  
  cat("## Remove Chimeras ---------------------------",  "\n")
  seqtab.nochim <- removeBimeraDenovo(seqtab.rightlen, method="consensus", 
                                      multithread=nCores, verbose = TRUE)
  cat("Count table no chimera dim:", dim(seqtab.nochim), "\n")
  # Fraction left with w/out chimera
  cat("Fraction without chimera:", sum(seqtab.nochim)/sum(seqtab.rightlen), "\n")
  saveRDS(seqtab.nochim, file.path(path.to.res, "seqtab_all_final.rds")) 
  
  cat("## Track reads after chimera check ---------------------------", "\n")
  track.file <- file.path(path.to.res, "track_all_pools.csv")
  track.df <- read.csv(file = track.file, row.names = 1)
  track <- cbind(track.df[rownames(seqtab.rightlen), ], 
                 rowSums(seqtab.rightlen), rowSums(seqtab.nochim))
  colnames(track) <- c("pool", "sample.id", "file.name", "input", "filtered", 
                       "merged", "tabled.rightlen", "nonchim")
  rownames(track) <- rownames(seqtab.rightlen)
  track.file <- file.path(path.to.res, "track_all_pools.csv")
  write.csv(track, file = track.file)
}

if(ASSIGN_TAXONOMY) {
  seqtab.nochim <- readRDS(file.path(path.to.res, "seqtab_all_final.rds"))
  cat("## Assign Taxonomy ---------------------------", "\n")
  # Use RDP and Silva reference databases
  refTaxa <- c(file.path(path.to.reference, "rdp_train_set_16.fa.gz"),
               file.path(path.to.reference, "silva_nr_v128_train_set.fa.gz"))
  
  refSpecies <- c(file.path(path.to.reference, "rdp_species_assignment_16.fa.gz"),
                  file.path(path.to.reference, "silva_species_assignment_v128.fa.gz"))
  names(refTaxa) <- names(refSpecies) <- c("RDP", "Silva")
  
  registerDoParallel(length(refTaxa))
  foreach(ref = names(refTaxa)) %dopar% {
    cat("## Running for ref:", ref,  "\n")
    taxa <- assignTaxonomy(seqs = colnames(seqtab.nochim), verbose = TRUE,
                           refFasta = refTaxa[[ref]], 
                           multithread = floor(nCores/length(refTaxa)))
    taxaSpecies <- addSpecies(taxa, refSpecies[[ref]], verbose=TRUE)
    ps <- phyloseq(otu_table(t(seqtab.nochim), taxa_are_rows = TRUE), 
                   tax_table(taxaSpecies))
    saveRDS(ps, file = file.path(path.to.res, paste0("phyloseq_", ref, ".rds")))
    
    taxaMultiSpecies <- addSpecies(taxa, refSpecies[[ref]], verbose=TRUE, 
                                   allowMultiple = TRUE)
    ps <- phyloseq(otu_table(t(seqtab.nochim), taxa_are_rows = TRUE), 
                   tax_table(taxaMultiSpecies))
    saveRDS(ps, file = file.path(path.to.res, paste0("phyloseq_", ref, 
                                                     "_multispecies.rds")))
  }
}


