library(testthat)
library(adegenet)
library(poppr)
library(openxlsx)
library(tidyr)

test_that("basic.qc runs without errors and returns expected output", {
   # Create dummy dataset
   dummy_data <- data.frame(
      Ind = c("Ind1", "Ind2", "Ind3", "Ind4"),
      Pop = c("Pop1", "Pop1", "Pop2", "Pop2"),
      SNP1 = c("A/A", "A/T", "T/T", "A/A"),
      SNP2 = c("G/G", "G/T", "T/T", "G/G")
   )
   write.csv(dummy_data, "test_input.csv", row.names = FALSE)
   
   # Run function on test data
   results <- suppressMessages(suppressWarnings(basic.qc("test_input.csv", remove.missing.samples = TRUE)))
   
   # Ensure exported file exists
   expect_true(file.exists("quality-control-results.xlsx"))
   
   # Load results from the saved file
   exported_data <- openxlsx::read.xlsx("quality-control-results.xlsx")
   expect_named(exported_data, names(results))
   
   # Clean up test files
   unlink("test_input.csv")
   unlink("quality-control-results.xlsx")
})