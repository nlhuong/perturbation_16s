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

#' Extract Function IDs across Species
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
#' Sys.setenv("MIDAS_DB" = "/scratch/users/kriss1/applications/MIDAS/database/midas_db_v1.2")
#' function_annotation(species_ids = "Acanthamoeba_endosymbiont_62344")
function_annotation <- function(species_ids = NULL, MIDAS_DB = NULL) {
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

#' Extract Interpretations for Function IDs
#'
#' This is just a convenience wrapper function for accessing the functional
#' information in the ontology subdirectory of the MIDAS reference database.
#'
#' @param function_ids [character vector] The function IDs to get some more
#'   human interpretable explanations for. For example, we want to know that
#'   FIG00000001 corresponds to Cysteine desulfurase.
#' @param MIDAS_DB [character] The directory containing the MIDAS database.
#'   Defaults to the environmental variable for the MIDAS database.
#' @examples
#' Sys.setenv("MIDAS_DB" = "/scratch/users/kriss1/applications/MIDAS/database/midas_db_v1.2")
#' function_explanation(c("FIG00000001", "FIG00000845"))
function_explanation <- function(function_ids, MIDAS_DB = NULL) {
  if (is.null(MIDAS_DB)) {
    MIDAS_DB <- Sys.getenv("MIDAS_DB")
  }

  ## paths to the data containing explanations
  ontology_dirs <- c(
    "ec" = file.path(MIDAS_DB, "ontologies", "ec.txt"),
    "figfam" = file.path(MIDAS_DB, "ontologies", "figfam.txt"),
    "go" = file.path(MIDAS_DB, "ontologies", "go.txt"),
    "path" = file.path(MIDAS_DB, "ontologies", "path.txt")
  )

  ## loop over ontology types and extract relevant functional info
  explanations <- list()
  for (i in seq_along(ontology_dirs)) {
    message("Reading from ", ontology_dirs[i])
    explanations[[i]] <- read_tsv(
      ontology_dirs[i],
      col_names = c("id", "function")
    ) %>%
      filter(id %in% function_ids) %>%
      mutate(ontology = basename(ontology_dirs[i]))
  }

  bind_rows(explanations)
}
