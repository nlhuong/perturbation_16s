#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Some functions to provide ontological (GO, KEGG, FIGFAM) annotation of the
## raw gene IDs output by MIDAS.
##
## author: sankaran.kris@gmail.com
## date: 12/31/2017

library("readr")
library("dplyr")

#' Extract Function Annotation DataFrame
#'
#' Given a list of species, extract the mapping between gene IDs and figfam / GO
#' / KEGG function annotations. This information is contained in the pan_genome
#' subdirectory of the MIDAS_DB, in the files with the name
#' centroid_functions.txt.gz (which are contained in each of the species
#' subdirectories).
#'
#' @param MIDAS_DB [character] The directory containing the MIDAS database.
#'   Defaults to the environmental variable for the MIDAS database.
#' @param species_ids [character vector] The list of species to focus on
#'   extracting centroid information about. A vector of species names to focus
#'   on extracting the functions for.
#' @return annotation [data.frame]
#' @examples
#' function_annotation(species_ids = "Acanthamoeba_endosymbiont_62344")
function_annotation <- function(MIDAS_DB = NULL, species_ids = NULL) {
  if (is.null(MIDAS_DB)) {
    MIDAS_DB <- Sys.getenv("MIDAS_DB")
  }

  ## get subdirectories of interest
  species_dirs <- list.files(file.path(MIDAS_DB, "pan_genomes"), full.name = TRUE)
  if (!is.null(species_dirs)) {
    species_dirs <- species_dirs[basename(species_dirs) %in% species_ids]
  }

  ## extract the functional information
  annotation <- list()
  for (i in seq_along(species_dirs)) {
    cat(sprintf(
      "Extracting functions for %s (%s / %s) \n",
      species_dirs[i], i, length(species_dirs)
    ))

    annotation[[i]] <- read_tsv(file.path(species_dirs[i], "centroid_functions.txt.gz"))
    annotation[[i]]$species_id <- basename(species_dirs[i])
    }

  bind_rows(annotation)
}
