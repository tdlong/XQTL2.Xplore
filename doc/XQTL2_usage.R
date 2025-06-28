## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
# # To update the genes data (from GTF file)
# source("data-raw/create_genes_object.R")
# devtools::document()
# devtools::install()
# data(dm6.ncbiRefSeq.genes)
# 
# # To update the variants data (from VCF files)
# source("data-raw/create_variants_object.R")
# devtools::document()
# devtools::install()
# data(dm6.variants)
# 
# # To update example datasets
# source("data-raw/create_example_data.R")
# devtools::document()
# devtools::install()
# data(zinc_hanson_pseudoscan)
# data(zinc_hanson_means)

## -----------------------------------------------------------------------------
library(XQTL2.Xplore)

# Load reference data
data(dm6.ncbiRefSeq.genes)
data(dm6.variants)

# Load example datasets
data(zinc_hanson_pseudoscan)
data(zinc_hanson_means)

## ----eval=FALSE---------------------------------------------------------------
# # Load your own QTL scan results
# my_qtl_data <- read.table("my_pseudoscan.txt", header = TRUE)
# 
# # Load your own frequency data
# my_freq_data <- read.table("my_meansBySample.txt", header = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# # Using example data
# XQTL_Manhattan_5panel(zinc_hanson_pseudoscan)
# 
# # Using your own data
# XQTL_Manhattan_5panel(my_qtl_data)
# 
# # Single Manhattan plot
# XQTL_Manhattan(zinc_hanson_pseudoscan)

## ----eval=FALSE---------------------------------------------------------------
# # Using example data
# A1 = XQTL_region(zinc_hanson_pseudoscan, "chr3R", 13000000, 15000000, "Wald_log10p")
# A2 = XQTL_change_average(zinc_hanson_means, "chr3R", 13000000, 15000000)
# A3 = XQTL_genes(dm6.ncbiRefSeq.genes, "chr3R", 13000000, 15000000)
# 
# # Using your own data
# A1 = XQTL_region(my_qtl_data, "chr3R", 13000000, 15000000, "Wald_log10p")
# A2 = XQTL_change_average(my_freq_data, "chr3R", 13000000, 15000000)
# A3 = XQTL_genes(dm6.ncbiRefSeq.genes, "chr3R", 13000000, 15000000)

## ----eval=FALSE---------------------------------------------------------------
# # Using example data
# out = XQTL_zoom(zinc_hanson_pseudoscan, "chr3R", 13000000, 15000000, left_drop = 3, right_drop = 3)
# 
# # Using your own data
# out = XQTL_zoom(my_qtl_data, "chr3R", 13000000, 15000000, left_drop = 3, right_drop = 3)

## ----eval=FALSE---------------------------------------------------------------
# # Example: Customizing a Manhattan plot
# p = XQTL_Manhattan(zinc_hanson_pseudoscan)
# p + theme_minimal() +
#     labs(title = "My Custom Title")

## ----eval=FALSE---------------------------------------------------------------
# # First find a peak
# out = XQTL_zoom(zinc_hanson_pseudoscan, "chr3R", 13000000, 15000000)
# 
# # Then plot the region around the peak
# A1 = XQTL_region(zinc_hanson_pseudoscan, out$chr, out$start, out$stop, "Wald_log10p")
# A2 = XQTL_change_average(zinc_hanson_means, out$chr, out$start, out$stop)
# A3 = XQTL_genes(dm6.ncbiRefSeq.genes, out$chr, out$start, out$stop)

