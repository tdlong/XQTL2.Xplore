# Load required libraries
library(XQTL2.Xplore)
library(tidyverse)
library(patchwork)

# Load the data
df1 = as_tibble(read.table("../bitnbobs/JUICE.pseudoscan.txt"))
df2 = as_tibble(read.table("../bitnbobs/JUICE.meansBySample.txt"))

# Load the genes data
data(dm6.ncbiRefSeq.genes)
data(dm6.variants)

# Get zoomed region
out = XQTL_zoom(df1, "chr3R", 13000000, 15000000, 3, 3)
out$plot  # perhaps adjust left and right drops until you are happy

# Create the four plots
A1 = XQTL_region(df1, out$chr, out$start, out$stop, "Wald_log10p")
A2 = XQTL_change_average(df2, out$chr, out$start, out$stop)
A3 = XQTL_genes(dm6.ncbiRefSeq.genes, out$chr, out$start, out$stop)
A4 = XQTL_variants(dm6.variants, out$chr, out$start, out$stop, df2)

# Display the combined plot
print(A1/A3/A4/A2) 