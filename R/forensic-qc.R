#' @title Basic QC of genotype files with population data
#' @description
#' This function performs quality control of genotype files with population data. 
#' @param input (character) The file path to the genetic dataset. Accepts csv or xlsx files. The first column lists the samples, the second column indicates their population, and the rest of the columns their genotype for a marker.
#' @param remove.missing.samples (logical) If TRUE, this will remove samples with more than a certain percentage of missing genotypes. Default is TRUE.
#' @param samples.percentage (integer) The missingness threshold for samples for removal.
#' @param remove.missing.snps (logical) If TRUE, this will remove SNPs with more than a certain percentage of missing genotypes. Default is FALSE.
#' @param snps.percentage (integer) The missingness threshold for snps for removal. Suggested only if not analyzing a certain set of SNPs. High missingness for those are informative.
#' @param remove.duplicates (logical) will remove samples with identical genotypes if set to TRUE. Default is TRUE.
#' @example quality.control("test_data/test.csv", remove.missing.samples = TRUE, samples.percentage = "80", remove.missing.snps = FALSE, remove.duplicates = TRUE)
#' @import ade4
#' @import tools
#' @import readxl
#' @import readr
#' @import dplyr
#' @import adegenet
#' @import poppr
#' @import tibble
#' @import factoextra
#' @import pegas
#' @import ape
#' @import openxlsx
#' @export

load_input_file <- function(input) {
   if (tools::file_ext(input) == "csv") {
      file <- readr::read_csv(input)
   } else if (tools::file_ext(input) == "xlsx") {
      file <- readxl::read_excel(input)
   } else {
      stop("Input file should be in csv or xlsx format.")
   }
   return(file)
}

clean_input_data <- function(file) {
   file <- lapply(file, function(x) gsub("|", "/", x, fixed = TRUE))
   file <- as.data.frame(file)
   
   file[is.na(file)] <- "N"
   file <- file %>%
      mutate(across(everything(), ~ case_when(
         . == "N/A" ~ "N", 
         . == "NA" ~ "N",
         TRUE ~ .x))) %>%
      rename(Ind = 1, Pop = 2)
   
   return(file)
}

convert_to_genind <- function(file) {
   ind <- as.character(file$Ind)
   pop <- as.character(file$Pop)
   fsnps_geno <- file[, 3:ncol(file)]
   
   fsnps_gen <- adegenet::df2genind(fsnps_geno, ind.names = ind, pop = pop, sep = "/", NA.char = "N", ploidy = 2, type = "codom")
   fsnps_gen@pop <- as.factor(file$Pop)
   
   return(fsnps_gen)
}

filter_missing_samples <- function(fsnps_gen, samples.percentage) {
   fsnps_gen_indmiss <- propTyped(fsnps_gen, by = "ind")
   missing_samples <- fsnps_gen_indmiss[which(fsnps_gen_indmiss < samples.percentage)]
   
   fsnps_gen_filtered <- poppr::missingno(fsnps_gen, type = "geno", cutoff = samples.percentage)
   return(list(filtered_data = fsnps_gen_filtered, missing_samples = missing_samples))
}

filter_missing_snps <- function(fsnps_gen, snps.percentage) {
   fsnps_gen_locmiss <- propTyped(fsnps_gen, by = "loc")
   missing_loci <- fsnps_gen_locmiss[which(fsnps_gen_locmiss < snps.percentage)]
   
   fsnps_gen_filtered <- poppr::missingno(fsnps_gen, type = "loci", cutoff = snps.percentage)
   return(list(filtered_data = fsnps_gen_filtered, missing_loci = missing_loci))
}

remove_duplicates <- function(fsnps_gen, file) {
   fsnps_gen_dups <- mlg.id(fsnps_gen)
   
   duplicates <- fsnps_gen_dups[lengths(fsnps_gen_dups) == 2]
   
   if (length(duplicates) == 0) {
      return(list(final_dataset = tibble::rownames_to_column(as.data.frame(fsnps_gen), "Ind"), duplicates = NULL))
   } else {
      dup2 <- data.frame(t(as.data.frame(duplicates)))
      dup2 <- dup2 %>% rename(Ind = 1, Ind2 = 2)
      
      file$Ind <- as.character(file$Ind)
      file_corr <- file %>% anti_join(dup2, by = "Ind")
      
      return(list(final_dataset = file_corr, duplicates = dup2))
   }
}

basic.qc <- function(input, remove.missing.samples = TRUE, samples.percentage = 0.80, remove.missing.snps = FALSE, snps.percentage = 0.90, remove.duplicates = TRUE) {
   if (!require("pacman")) install.packages("pacman")
   pacman::p_load(ade4, adegenet, factoextra, pegas, poppr, ape, dplyr, readr, openxlsx)
   
   # Load and clean data
   file <- load_input_file(input)
   file <- clean_input_data(file)
   
   # Convert to genind format
   fsnps_gen <- convert_to_genind(file)
   
   # Identify and filter missing samples
   if (remove.missing.samples) {
      sample_results <- filter_missing_samples(fsnps_gen, samples.percentage)
      fsnps_gen <- sample_results$filtered_data
      missing_samples <- sample_results$missing_samples
   } else {
      missing_samples <- NULL
   }
   
   # Identify and filter missing SNPs
   if (remove.missing.snps) {
      snps_results <- filter_missing_snps(fsnps_gen, snps.percentage)
      fsnps_gen <- snps_results$filtered_data
      missing_loci <- snps_results$missing_loci
   } else {
      missing_loci <- NULL
   }
   
   # Identify and remove duplicates
   if (remove.duplicates) {
      dup_results <- remove_duplicates(fsnps_gen, file)
      final_dataset <- dup_results$final_dataset
      duplicates <- dup_results$duplicates
   } else {
      final_dataset <- tibble::rownames_to_column(as.data.frame(fsnps_gen), "Ind")
      duplicates <- NULL
   }
   
   # Export results
   datasets <- list(
      "QC Dataset" = final_dataset,
      "Duplicates" = duplicates,
      "Samples with High Missingness" = missing_samples,
      "Loci with High Missingness" = missing_loci
   )
   
   openxlsx::write.xlsx(datasets, file = "quality-control-results.xlsx")
}