## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 8,
  fig.height = 6
)

## ----load_data----------------------------------------------------------------
library(XQTL2.Xplore)

# Load example datasets
data(zinc_hanson_pseudoscan)  # QTL scan results
data(zinc_hanson_means)       # Founder frequency data

# Load reference data for genes and variants
data(dm6.ncbiRefSeq.genes)    # Gene annotations
data(dm6.variants)           # Variant data

# Display basic information about the datasets
cat("ZINC Hanson pseudoscan data dimensions:", dim(zinc_hanson_pseudoscan), "\n")
cat("ZINC Hanson means data dimensions:", dim(zinc_hanson_means), "\n")
cat("Available chromosomes:", unique(zinc_hanson_pseudoscan$chr), "\n")

## ----manhattan_plots----------------------------------------------------------
# 5-panel Manhattan plot with physical distance (Mb) - useful for precise peak localization
XQTL_Manhattan_5panel(zinc_hanson_pseudoscan, cM = FALSE)

# 5-panel Manhattan plot with genetic distance (cM) - useful for broad peak identification
XQTL_Manhattan_5panel(zinc_hanson_pseudoscan, cM = TRUE)

# Traditional Manhattan plot with physical distance
XQTL_Manhattan(zinc_hanson_pseudoscan, cM = FALSE, color_scheme = "UCI")

# Traditional Manhattan plot with genetic distance
XQTL_Manhattan(zinc_hanson_pseudoscan, cM = TRUE)

## ----frequency_analysis-------------------------------------------------------
# Basic frequency change analysis
XQTL_change_average(zinc_hanson_means, "chr3R", 18000000, 20000000)

# Frequency change with reference strain highlighting
XQTL_change_average(zinc_hanson_means, "chr3R", 18000000, 20000000, reference_strain = "A1")

# Frequency changes by replicate
XQTL_change_byRep(zinc_hanson_means, "chr3R", 18000000, 20000000)

## ----peak_refinement----------------------------------------------------------
# Find and refine peak boundaries
# Arguments: data, chromosome, start, stop, left_drop, right_drop
out <- XQTL_zoom(zinc_hanson_pseudoscan, "chr3R", 18000000, 20000000, 3, 3)

# View the refined peak
out$plot

# Display the new interval boundaries
cat("Refined interval:", out$chr, ":", out$start, "-", out$stop, "\n")

## ----regional_analysis--------------------------------------------------------
# Create individual plots for the refined region
A1 <- XQTL_region(zinc_hanson_pseudoscan, out$chr, out$start, out$stop, "Wald_log10p")
A2 <- XQTL_change_average(zinc_hanson_means, out$chr, out$start, out$stop)
A2b <- XQTL_change_average(zinc_hanson_means, out$chr, out$start, out$stop, plotSelection = TRUE)
A3 <- XQTL_genes(dm6.ncbiRefSeq.genes, out$chr, out$start, out$stop)
A4 <- XQTL_variantsByFounder(dm6.variants, out$chr, out$start, out$stop, zinc_hanson_means)
A5 <- XQTL_SVBySize(dm6.variants, out$chr, out$start, out$stop, zinc_hanson_means)

# Combine plots using patchwork
library(patchwork)
A1 / A3 / A2
A4 / A5

## ----multi_panel--------------------------------------------------------------
# Create publication-ready 5-panel plot
XQTL_5panel_plot(zinc_hanson_pseudoscan, zinc_hanson_means, 
                 dm6.variants, dm6.ncbiRefSeq.genes, 
                 out$chr, out$start, out$stop)

## ----second_peak--------------------------------------------------------------
# Explore a peak on chromosome 2L
out2 <- XQTL_zoom(zinc_hanson_pseudoscan, "chr2L", 6000000, 7500000, 2, 2)
out2$plot

# Create multi-panel plot for this region
myplot <- XQTL_5panel_plot(zinc_hanson_pseudoscan, zinc_hanson_means, 
                           dm6.variants, dm6.ncbiRefSeq.genes, 
                           out2$chr, out2$start, out2$stop)

# Save the plot (uncomment to save)
# png("combined_plot_test.png", height = 10, width = 7, units = "in", res = 300)
# print(myplot)
# dev.off()

