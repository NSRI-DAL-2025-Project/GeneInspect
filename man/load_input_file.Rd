% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/explore.pca.R, R/forensic-qc.R, R/pop.stat.R
\name{load_input_file}
\alias{load_input_file}
\title{Exploratory PCA}
\usage{
load_input_file(input)

load_input_file(input)

load_input_file(input)
}
\arguments{
\item{input}{(character) The file path to the genetic dataset. Accepts .csv or .xlsx files.}

\item{default.colors.labels}{(logical) If TRUE, default colors from Set1 of ColorBrewer will be used. Default is TRUE.}

\item{pca.labels}{(character) The file path to the text file that should be supplied if there is an order of the labels preferred. Required if default.colors.labels is FALSE.}

\item{color.palette}{(character) The file path to the text file containing hex codes that should be supplied if the labels will be color coded.}

\item{set.size}{(logical) If TRUE, user will set the size of the PNG output. Default is FALSE.}

\item{width}{(integer) the width size of the PNG file. Should be indicated if set.size is TRUE.}

\item{height}{(integer) the height of the PNG file. Should be indicated if set.size is TRUE.}

\item{add.pc}{(logical) If TRUE, additional principal components will be plotted. Default is FALSE.}

\item{add.pc.x}{(character) supply the PC in the X axis. It should follow the format "PC1", "PC2", "PC3", ..., "PC6".}

\item{add.pc.y}{(character) supply the PC in the Y axis. It should follow the format "PC1", "PC2", "PC3", ..., "PC6".}

\item{remove.missing.samples}{(logical) If TRUE, this will remove samples with more than a certain percentage of missing genotypes. Default is TRUE.}

\item{samples.percentage}{(integer) The missingness threshold for samples for removal.}

\item{remove.missing.snps}{(logical) If TRUE, this will remove SNPs with more than a certain percentage of missing genotypes. Default is FALSE.}

\item{snps.percentage}{(integer) The missingness threshold for snps for removal. Suggested only if not analyzing a certain set of SNPs. High missingness for those are informative.}

\item{remove.duplicates}{(logical) will remove samples with identical genotypes if set to TRUE. Default is TRUE.}

\item{markers}{(integer) (Optional) Number of genetic markers to analyze. Defaults to all available markers.}

\item{order.label}{(logical) If TRUE, population labels are ordered in custom. Default is FALSE.}

\item{pop.labels}{(character) (Optional) Custom population labels for plotting and grouping.}
}
\description{
Function for the generation of PCA plots

This function performs quality control of genotype files with population data.

pop.stats() is an R function designed for population genetic analysis. It takes in genetic data as input and processes basic and keys statistical measures and exports the results.
}
