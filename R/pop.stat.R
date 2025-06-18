#' @title pop.stats Documentation
#' @description pop.stats() is an R function designed for population genetic analysis. It takes in genetic data as input and processes basic and keys statistical measures and exports the results.
#' @param input (character) The file path to the genetic dataset. Accepts .csv or .xlsx files.
#' @param markers (integer) (Optional) Number of genetic markers to analyze. Defaults to all available markers.
#' @param order.label (logical) If TRUE, population labels are ordered in custom. Default is FALSE.
#' @param pop.labels (character) (Optional) Custom population labels for plotting and grouping.
#' @example pop.stats("test_data/test.csv", markers = 56, order.labels = FALSE)
#' @import ade4
#' @import adegenet
#' @import factoextra
#' @import pegas
#' @import poppr
#' @import ape
#' @import dplyr
#' @import hierfstat
#' @import reshape2
#' @import ggplot2
#' @import ggrepel
#' @import ggpubr
#' @import RColorBrewer
#' @import scales
#' @import devtools
#' @import miscTools
#' @import stringr
#' @import lattice
#' @import vegan
#' @import readr
#' @import openxlsx
#' @export
load_input_file <- function(input) {
   if (tools::file_ext(input) == "csv") {
      return(readr::read_csv(input))
   } else if (tools::file_ext(input) == "xlsx") {
      return(readxl::read_excel(input))
   } else {
      stop("Input file should be in csv or xlsx format.")
   }
}

clean_input_data <- function(file) {
   file <- lapply(file, function(x) gsub("|", "/", x, fixed = TRUE))
   file <- as.data.frame(file)
   file[is.na(file)] <- "N"
   
   file <- file %>%
      mutate(across(everything(), ~ case_when(. == "N/A" ~ "N", . == "NA" ~ "N", TRUE ~ .x))) %>%
      rename(Ind = 1, Pop = 2)
   
   file <- as.data.frame(file)
   
   return(file)
}

convert_to_genind <- function(file) {
   ind <- as.character(file$Ind)
   pop <- as.character(file$Pop)
   fsnps_geno <- file[, 3:ncol(file)]
   
   fsnps_gen <- adegenet::df2genind(fsnps_geno, 
                                    ind.names = ind, 
                                    pop = pop, 
                                    sep = "/", 
                                    NA.char = "N", 
                                    ploidy = 2, 
                                    type = "codom")
   
   fsnps_gen@pop <- as.factor(file$Pop)
   
   return(fsnps_gen)
}

# Need to correct, genind2hierfstat is outputting an error 
compute_population_stats <- function(fsnps_gen) {
   #fsnps_gen@pop <- as.factor(fsnps_gen@pop)  # Ensure populations are correctly stored
   
   # Private Alleles (Matrix/List for export)
   priv_alleles <- poppr::private_alleles(fsnps_gen)
   if (is.null(priv_alleles)) priv_alleles <- list(message = "No private alleles detected")  # Ensure it's a list
   
   # Allelic Richness (Matrix/List for export)
   mar_matrix <- allelic.richness(genind2hierfstat(fsnps_gen))$Ar %>%
      apply(MARGIN = 2, FUN = mean) %>% 
      round(digits = 3)
   mar_list <- as.list(mar_matrix)

   # Heterozygosities (Data Frame for ggplot)
   basic_fsnps <- hierfstat::basic.stats(fsnps_gen, diploid = TRUE)
   Ho_fsnps <- apply(basic_fsnps$Ho, 2, mean, na.rm = TRUE) %>% round(2)
   He_fsnps <- apply(basic_fsnps$Hs, 2, mean, na.rm = TRUE) %>% round(2)
   
   Het_fsnps_df <- data.frame(Pop = names(Ho_fsnps), Ho = Ho_fsnps, He = He_fsnps) %>%
      tidyr::pivot_longer(cols = c("Ho", "He"), names_to = "Variable", values_to = "Value")
   
   return(list(stats_matrix = list(priv_alleles = priv_alleles, mar_list = mar_list), 
               stats_df = Het_fsnps_df))
   
}

compute_hardy_weinberg <- function(fsnps_gen) {
   # Hardy-Weinberg Equilibrium (List for export)
   fsnps_hwe <- as.numeric(round(pegas::hw.test(fsnps_gen, B = 1000), 6)) 
   
   # Chi-square test (Matrix for export, Data Frame for ggplot)
   fsnps_hwe_test <- data.frame(sapply(adegenet::seppop(fsnps_gen), 
                                       function(ls) pegas::hw.test(ls, B=0)[,3]))
   
   fsnps_hwe_chisq_matrix <- as.matrix(fsnps_hwe_test)  # Convert to matrix
   fsnps_hwe_chisq_df <- as.data.frame(t(fsnps_hwe_chisq_matrix)) %>% tibble::rownames_to_column("Population")
   
   return(list(hw_matrix = list(hwe = fsnps_hwe, hwe_chisq_matrix = fsnps_hwe_chisq_matrix), 
               hw_df = fsnps_hwe_chisq_df))
}

compute_fst <- function(fsnps_gen) {
   # Compute Fst values as a matrix
   fsnps_fst_matrix <- as.matrix(hierfstat::genet.dist(fsnps_gen, method = "WC84") %>% round(3))
   fst_list <- as.list(fsnps_fst_matrix)  # Convert matrix to list
   if (length(fst_list) == 0) fst_list <- list(message = "No Fst values calculated")  # Ensure it's never NULL
   
   # Convert matrix into a named list
   fst_list <- as.list(fsnps_fst_matrix)  # Ensure it's list format
   
   # Convert matrix into a tidy data frame for ggplot
   fst_df <- as.data.frame(fsnps_fst_matrix) %>%
      tibble::rownames_to_column(var = "Site1") %>% 
      tidyr::pivot_longer(cols = -Site1, names_to = "Site2", values_to = "Fst")
   
   return(list(fst_matrix = fst_list, fst_df = fst_df))  # Returns correctly formatted values
}


plot_heterozygosity <- function(Het_fsnps_df) {
   ggplot(Het_fsnps_df, aes(x = Pop, y = Value, fill = Variable)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.6), colour = "black") +
      scale_y_continuous(expand = c(0,0), limits = c(0,0.50)) +
      scale_fill_manual(values = c("royalblue", "#bdbdbd"), labels = c(expression(italic("H")[o]), expression(italic("H")[e]))) +
      labs(y = "Heterozygosity") +
      theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, face = "bold"))
   
   ggsave(filename = "heterozygosity_plot.tiff", width = 9, dpi = 300)
}

plot_fst_heatmap <- function(fst_df) {
   ggplot(fst_df, aes(x = Site1, y = Site2, fill = Fst)) +
      geom_tile(colour = "black") +
      geom_text(aes(label = Fst), color="black", size = 3) +
      scale_fill_gradient2(low = "blue", mid = "pink", high = "red", midpoint = max(fst_df$Fst) / 2) +
      theme(axis.text = element_text(colour = "black", size = 8, face = "bold"))
   
   ggsave(filename = "fst_plot.tiff", width = 8, height = 6, dpi = 600)
}


export_results <- function(stats_matrix, hw_matrix, fst_matrix) {
   # Verify fst_matrix is a proper list
   if (!is.list(fst_matrix)) {
      fst_matrix <- as.list(fst_matrix)  # Convert if necessary
   }
   
   priv_alleles <- as.data.frame(stats_matrix$priv_alleles)
   mar <- as.data.frame(stats_matrix$mar_list)
   heterozygosities <- as.data.frame(stats_matrix$stats_df)
   hwe <- as.data.frame(hw_matrix$hwe)
   chisquare <- as.data.frame(hw_matrix$hwe_chisq_matrix)
   fst <- as.data.frame(fst_matrix$fst_matrix)
   
   datasets <- list(
      "Private Alleles" = priv_alleles,
      "Mean allelic richness" = mar,
      "Heterozygosities" = heterozygosities,
      "Hardy-Weinberg Equilibrium" = hwe,
      "Chi-square test" = chisquare,
      "Fst values" = fst
   )
   
   openxlsx::write.xlsx(datasets, file = "population-statistics-results.xlsx")
   
   # Validate file contents in the test
   exported_data <- openxlsx::read.xlsx("quality-control-results.xlsx")
   expect_named(exported_data, names(datasets))  # Ensure correct labels
}

pop.stats <- function(input, markers = NULL, order.labels = FALSE, pop.labels = NULL) {
   if (!require("pacman")) install.packages("pacman")
   pacman::p_load(ade4, adegenet, factoextra, pegas, poppr, ape, dplyr, hierfstat, reshape2, ggplot2, ggrepel, ggpubr,
                  RColorBrewer, scales, devtools, miscTools, stringr, lattice, vegan, readr, openxlsx)
   
   file <- load_input_file(input)
   file <- clean_input_data(file)
   fsnps_gen <- convert_to_genind(file)
   
   stats_results <- compute_population_stats(fsnps_gen)
   hw_results <- compute_hardy_weinberg(fsnps_gen)
   fst_results <- compute_fst(fsnps_gen)
   
   # Export statistical results
   export_results(stats_results$stats_matrix, hw_results$hw_matrix, fst_results$fst_matrix)
   
   # Visualize heterozygosity & Fst heatmap
   plot_heterozygosity(stats_results$stats_df)
   plot_fst_heatmap(fst_results$fst_df)
   
   return(list(stats_results = stats_results$stats_matrix, hw_results = hw_results$hw_matrix, fst_results = fst_results$fst_matrix))
}


