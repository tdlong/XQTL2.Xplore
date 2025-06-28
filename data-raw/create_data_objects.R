# This script creates the data objects that will be included in the package
# Run this script from the package root directory
# 
# This is a wrapper that calls the two working scripts:
# 1. create_genes_object.R - processes GTF file (fast)
# 2. create_variants_object.R - processes VCF files (fast)

message("Creating data objects using working scripts...")

# Create the genes object using the working script
message("Step 1: Creating genes object...")
source("data-raw/create_genes_object.R")
message("Genes object created successfully!")

# Create the variants object using the working script  
message("Step 2: Creating variants object...")
source("data-raw/create_variants_object.R")
message("Variants object created successfully!")

message("All data objects created successfully!")
message("Objects created:")
message("- dm6.ncbiRefSeq.genes")
message("- dm6.variants") 