library(testthat)
library(adegenet)
library(poppr)
library(openxlsx)
library(tidyr)

test_that("pop.stats runs without errors and returns expected output", {
   # Create dummy dataset
   dummy_data <- data.frame(
      Ind = c("Ind1", "Ind2", "Ind3", "Ind4"),
      Pop = c("Pop1", "Pop1", "Pop2", "Pop2"),
      SNP1 = c("A/A", "A/T", "T/T", "A/A"),
      SNP2 = c("G/G", "G/T", "T/T", "G/G")
   )
   write.csv(dummy_data, "test_input.csv", row.names = FALSE)
   
   # Run function on test data
   results <- suppressMessages(suppressWarnings(pop.stats("test_input.csv", markers = 2)))  
   
   # Validate structure
   expect_type(results, "list")
   expect_named(results, c("stats_results", "hw_results", "fst_results"))
   
   # Validate data types
   expect_type(results$stats_results$priv_alleles, "list")  
   expect_type(results$stats_results$mar_list, "list")  
   
   expect_type(results$hw_results$hwe, "numeric")  
   expect_type(results$hw_results$hwe_chisq_matrix, "matrix")  
   
   expect_type(results$fst_results$fst_matrix, "list")  
   expect_true(length(results$fst_results$fst_matrix) > 0)  
   
   # Validate exported file
   expect_true(file.exists("quality-control-results.xlsx"))
   exported_data <- openxlsx::read.xlsx("quality-control-results.xlsx")
   expect_named(exported_data, c("Private Alleles", "Mean allelic richness", "Hardy-Weinberg Equilibrium", "Chi-square test", "Fst values"))
   
   # Validate ggplot-compatible outputs
   expect_s3_class(results$stats_results$stats_df, "data.frame")  
   expect_s3_class(results$fst_results$fst_df, "data.frame")  
   
   # Clean up
   unlink("test_input.csv")
   unlink("quality-control-results.xlsx")
})