library(testthat)
library(adegenet)
library(RColorBrewer)

# Test: Load Input File
test_that("load_input_file loads valid files", {
   csv_file <- system.file("test_data", "test.csv", package = "basics.forensics")
   
   expect_false(csv_file == "", "Test file not found.")  # Check if file exists
   expect_s3_class(load_input_file(csv_file), "data.frame")
})

# Test: Compute PCA
test_that("compute_pca returns correct data structure", {
   dummy_genind <- df2genind(data.frame(SNP1 = c("A/T", "T/T")), 
                             ind.names = c("Ind1", "Ind2"), 
                             pop = c("Pop1", "Pop2"), 
                             sep = "/", ploidy = 2, type = "codom")
   
   pca_results <- compute_pca(dummy_genind)
   
   expect_type(pca_results, "list")  # Ensure correct return type
   expect_named(pca_results, c("pca1", "percent", "ind_coords", "centroid"))  # Validate structure
   
   expect_true(is.numeric(pca_results$percent), "Percent values should be numeric")
   expect_true(is.data.frame(pca_results$ind_coords), "Individual coordinates must be a dataframe")
   expect_true(is.data.frame(pca_results$centroid), "Centroid data must be a dataframe")
})

# Test: Color and Labels Handling
test_that("get_colors_labels produces correct list format", {
   dummy_genind <- df2genind(data.frame(SNP1 = c("A/T", "T/T")), 
                             ind.names = c("Ind1", "Ind2"), 
                             pop = c("Pop1", "Pop2"), 
                             sep = "/", ploidy = 2, type = "codom")
   
   labels_colors <- get_colors_labels(dummy_genind, default.colors.labels = TRUE)
   
   expect_type(labels_colors, "list")  # Ensure it's a list
   expect_named(labels_colors, c("labels", "colors"))  # Validate expected names
})

# Test: PCA Plot Generation
test_that("plot_pca generates expected plots", {
   dummy_data <- data.frame(PC1 = rnorm(10), PC2 = rnorm(10), Site = rep(c("Pop1", "Pop2"), each = 5))
   dummy_centroid <- aggregate(cbind(PC1, PC2) ~ Site, data = dummy_data, FUN = mean)
   dummy_percent <- rep(10, 7)  # Simulated percentage variance
   dummy_colors <- list(labels = c("Pop1", "Pop2"), colors = c("#E41A1C", "#377EB8"))
   
   expect_type(dummy_colors, "list")  # Ensure dummy_colors is structured correctly
   
   expect_silent(plot_pca(dummy_data, dummy_centroid, dummy_percent, dummy_colors, filename = "test_plot.png"))
   expect_true(file.exists("test_plot.png"), "Plot file should be created")
   
})