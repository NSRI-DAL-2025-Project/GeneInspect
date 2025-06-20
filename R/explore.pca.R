#' @title Exploratory PCA
#' @description Function for the generation of PCA plots
#' @param input (character) The file path to the genetic dataset. Accepts csv or xlsx files. Genotypes in columns and samples in rows. First column should contain the samples, second column the population, and the rest of the columns are the markers.
#' @param default.colors.labels (logical) If TRUE, default colors from Set1 of ColorBrewer will be used. Default is TRUE.
#' @param pca.labels (character) The file path to the text file that should be supplied if there is an order of the labels preferred. Required if default.colors.labels is FALSE.
#' @param color.palette (character) The file path to the text file containing hex codes that should be supplied if the labels will be color coded.
#' @param set.size (logical) If TRUE, user will set the size of the PNG output. Default is FALSE.
#' @param width (integer) the width size of the PNG file. Should be indicated if set.size is TRUE.
#' @param height (integer) the height of the PNG file. Should be indicated if set.size is TRUE.
#' @param add.pc (logical) If TRUE, additional principal components will be plotted. Default is FALSE.
#' @param add.pc.x (character) supply the PC in the X axis. It should follow the format "PC1", "PC2", "PC3", ..., "PC6".
#' @param add.pc.y (character) supply the PC in the Y axis. It should follow the format "PC1", "PC2", "PC3", ..., "PC6".
#' @example explore.pca("test_data/test.csv", default.colors.labels = FALSE, set.size = FALSE, add.pc = TRUE, add.pc.x = "PC5", add.pc.y = "PC6")
#' @example explore.pca("test_data/test.csv", default.colors.labels = "labels.txt", color.palette = "colors.txt", set.size = TRUE, width = 14, height = 12, add.pc = FALSE)
#' @import tools
#' @import dplyr
#' @import adegenet
#' @import ggplot2
#' @import RColorBrewer
#' @export
#' 
# Load the input file (CSV or XLSX)
load_input_file <- function(input) {
   if (tools::file_ext(input) == "csv") {
      return(readr::read_csv(input))
   } else if (tools::file_ext(input) == "xlsx") {
      return(readxl::read_excel(input))
   } else {
      stop("Input file should be in csv or xlsx format.")
   }
}

#' @export
# Clean the input data
clean_input_data <- function(file) {
   file <- lapply(file, function(x) gsub(pattern = "|", replacement = "/", x = x, fixed = TRUE))
   file <- as.data.frame(file)
   file[is.na(file)] <- "N"
   file <- file %>% mutate(across(everything(), ~ case_when(
      . == "N/A" ~ "N", 
      . == "NA" ~ "N", 
      . == "." ~ "N",
      TRUE ~ .x)))
   file <- file %>% rename(Ind = 1, Pop = 2)
   return(file)
}

#' @export
# Convert cleaned data to genind format
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

#' @export
# Perform PCA computation
compute_pca <- function(fsnps_gen) {
   x <- tab(fsnps_gen, NA.method = "mean")
   set.seed(9999)
   pca1 <- dudi.pca(x, scannf = FALSE, scale = FALSE, nf = min(7, ncol(x)))
   
   percent <- pca1$eig / sum(pca1$eig) * 100
   ind_coords <- as.data.frame(pca1$li)  # Convert matrix to data frame
   colnames(ind_coords) <- paste0("PC", seq_len(ncol(ind_coords)))
   
   ind_coords$Ind <- indNames(fsnps_gen)
   ind_coords$Site <- fsnps_gen@pop
   
   centroid <- aggregate(ind_coords[, -c(ncol(ind_coords), ncol(ind_coords)-1)], 
                         by = list(ind_coords$Site), 
                         FUN = mean)
   colnames(centroid)[1] <- "Site"
   centroid <- as.data.frame(centroid)  # Convert matrix to data frame
   
   return(list(pca1 = pca1, percent = percent, ind_coords = ind_coords, centroid = centroid))
}

#' @export
get_colors_labels <- function(fsnps_gen, default.colors.labels, pca.labels = NULL, color.palette = NULL) {
   if (default.colors.labels) {
      labels <- levels(as.factor(fsnps_gen@pop))
      colors <- brewer.pal(nPop(fsnps_gen), "Set1")
   } else {
      labels <- read.csv(pca.labels, header = FALSE)$V1
      colors <- read.csv(color.palette, header = FALSE)$V1
      
      if (length(labels) != length(colors)) {
         stop("The number of labels and colors must match.")
      }
   }
   
   return(list(labels = labels, colors = colors))  # Ensure named list
}

#' @export
# Plotting 
plot_pca <- function(ind_coords, centroid, percent, labels_colors, filename, width = 8, height = 8, pc_x = 1, pc_y = 2) {
   # Convert matrices to data frames if needed
   if (!is.data.frame(ind_coords)) ind_coords <- as.data.frame(ind_coords)
   if (!is.data.frame(centroid)) centroid <- as.data.frame(centroid)
   
   if (!is.list(labels_colors)) {
      stop("labels_colors must be a named list containing 'labels' and 'colors'.")
   }
   
   # Ensure Site column exists and matches colors
   ind_coords$Site <- factor(ind_coords$Site, levels = labels_colors$labels)
   colors_named <- setNames(labels_colors$colors, labels_colors$labels)
   
   # Labels for axes
   xlab <- paste("PC", pc_x, " (", format(round(percent[pc_x], 1), nsmall = 1), "%)", sep = "")
   ylab <- paste("PC", pc_y, " (", format(round(percent[pc_y], 1), nsmall = 1), "%)", sep = "")
   
   # Create plot
   plot <- ggplot(data = ind_coords, aes_string(x = paste0("PC", pc_x), y = paste0("PC", pc_y))) +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      geom_point(aes(fill = Site), shape = 21, size = 4, show.legend = FALSE) +
      geom_label_repel(data = centroid, aes_string(x = paste0("PC", pc_x), y = paste0("PC", pc_y), label = "Site"), 
                       size = 4, show.legend = FALSE, max.overlaps = Inf) +
      scale_fill_manual(values = colors_named) +
      scale_colour_manual(values = colors_named) +
      labs(x = xlab, y = ylab)
   
   # Save plot
   ggsave(filename = filename, plot = plot, width = width, height = height, dpi = 600)
}

#' @export
# Main function
explore.pca <- function(input, default.colors.labels = TRUE, pca.labels = NULL, color.palette = NULL, set.size = FALSE, width = NULL, height = NULL, add.pc = FALSE, add.pc.x = NULL, add.pc.y = NULL) {
   # Step 1: Load and clean data
   file <- load_input_file(input)
   file <- clean_input_data(file)
   
   # Step 2: Convert to genind format
   fsnps_gen <- convert_to_genind(file)
   
   # Step 3: Perform PCA analysis
   pca_results <- compute_pca(fsnps_gen)
   
   # Step 4: Determine colors and labels (either default or user-specified)
   labels_colors <- get_colors_labels(fsnps_gen, default.colors.labels, pca.labels, color.palette)
   
   # Step 5: Validate width and height if set.size = TRUE
   if (set.size) {
      if (is.null(width) || is.null(height)) {
         stop("If set.size is TRUE, both width and height must be provided.")
      }
   } else {
      width <- 8
      height <- 8
   }
   
   # Step 6: Generate PCA plot with default PCs
   plot_pca(pca_results$ind_coords, pca_results$centroid, pca_results$percent, labels_colors, filename = "pca.png", width = width, height = height, pc_x = 1, pc_y = 2)
   
   # Step 7: Additional PCA plot if user specifies PCs
   if (add.pc) {
      plot_pca(pca_results$ind_coords, pca_results$centroid, pca_results$percent, labels_colors, filename = "pca2.png", width = width, height = height, pc_x = add.pc.x, pc_y = add.pc.y)
   }
   
   print("PCA analysis and plotting completed!")
}