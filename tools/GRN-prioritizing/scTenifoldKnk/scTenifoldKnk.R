#!/usr/bin/env Rscript

library(scTenifoldKnk)
library(dplyr)
library(optparse)
library(tidyverse)

#
# Define the arguments passed to the scTenifoldKnk tool
#
option_list <- list(
  make_option(c("-M", "--countMatrix"), 
              type = "character", 
              default = "countMatrix.csv", 
              help = "Path to the input count/expression matrix as a CSV.", 
              metavar = "character"),
  make_option(c("-P", "--perturbedGenes"), 
              type = "character", 
              default = "perturbedGenes.csv", 
              help = "Path to the output perturbed genes file [default %default]", 
              metavar = "character"),
  make_option(c("-G", "--geneKnockout"), 
              type = "character", 
              default = NULL,
              help = "A single gene or a comma-separated list of genes to be knocked out. For example gene1,gene2,gene3."),
  make_option(c("-S", "--qcMinLSize"), 
              type = "integer", 
              default = 0,
              help = "An integer value. Defines the minimum library size required for a cell to be included in the analysis."),
  make_option(c("-s", "--seed"), 
              type = "integer", 
              default = 42,
              help = "Set the seed for reproducibility."),
  make_option(c("-v", "--verbose"), 
              action = "store_true", 
              default = TRUE,
              help = "Print extra output."),
)

#
# Parse the arguments
#
opt <- parse_args(OptionParser(option_list = option_list))

#
# Set the seed.
#
set.seed(opt$seed)

#
# Check whether the count matrix file exists
#
if (!is.null(opt$countMatrix)) {
  if (!file.exists(opt$countMatrix)) {
    stop("Error. The specified input file does not exist.")
  }
} else {
  stop("Error. No input file provided. Use --countMatrix to specify the file path.")
}


#
# Process the geneKnockout argument
#
if (!is.null(opt$geneKnockout)) {
  gKO <- strsplit(opt$geneKnockout, ",")[[1]]
  gKO <- trimws(gKO)
  
  cat("Genes to knock out:", paste(gKO, collapse = ", "), "\n")
} else {
  stop("Error. No genes specified for knockout. Use --geneKnockout to specify as gene or genes in comma separated manner.\n")
}

#
# Read the count/expression matrix from a CSV file 
# that came from the AnnData object and format it
#
message("Reading the count matrix...")
countMatrix <- read.csv(opt$countMatrix,
                        header = TRUE)
rownames(countMatrix) <- countMatrix[ , 1]
countMatrix[ , 1] <- NULL
countMatrix <- t(countMatrix) %>%
  as.data.frame()

#
# Run scTenifoldKnk
#
message("Running scTenifoldKnk...")
res <- scTenifoldKnk(countMatrix = countMatrix, 
                     gKO = gKO, 
                     qc_minLSize = opt$qcMinLSize)

#
# Extract the data for comparison with ground truth
#
perturbedGenes <- res$diffRegulation %>%
  arrange(desc(Z)) %>%
  slice(c(1:5, (n() - 4):n())) %>%
  select(gene)
  
#
# Save the perturbed gene CSV
#
message("Saving perturbed genes...")
write.csv(x = perturbedGenes,
          file = opt$perturbedGenes,
          row.names = FALSE)
