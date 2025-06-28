# Installation Guide for XQTL2.Xplore

## Prerequisites

Before installing XQTL2.Xplore, make sure you have the following R packages installed:

```r
# Install required packages
install.packages(c("devtools", "ggplot2", "dplyr", "tidyr", "patchwork", "tibble", "grid", "rlang"))

# Install Bioconductor packages (if not already installed)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("vcfR")
```

## Installation from GitHub

### Method 1: Using devtools (Recommended for RStudio)

```r
# Install from GitHub with vignettes (important for RStudio users!)
devtools::install_github("tdlong/XQTL2.Xplore", build_vignettes = TRUE)

# Load the package
library(XQTL2.Xplore)
```

**Note:** The `build_vignettes = TRUE` parameter is crucial for RStudio users to access the included tutorials.

### Method 2: Using remotes

```r
# Install remotes if not already installed
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

# Install from GitHub with vignettes
remotes::install_github("tdlong/XQTL2.Xplore", build_vignettes = TRUE)

# Load the package
library(XQTL2.Xplore)
```

### Method 3: Manual Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/tdlong/XQTL2.Xplore.git
   cd XQTL2.Xplore
   ```

2. Install in R:
   ```r
   # Set working directory to the package folder
   setwd("path/to/XQTL2.Xplore")
   
   # Install the package with vignettes
   devtools::install(build_vignettes = TRUE)
   
   # Load the package
   library(XQTL2.Xplore)
   ```

## Quick Start

After installation, you can start using the package immediately:

```r
# Load the package
library(XQTL2.Xplore)

# Load example datasets
data(zinc_hanson_pseudoscan)
data(zinc_hanson_means)

# Create a Manhattan plot
XQTL_Manhattan_5panel(zinc_hanson_pseudoscan, cM = FALSE)
```

## Vignettes

Access the included tutorials:

```r
# View available vignettes
browseVignettes("XQTL2.Xplore")

# Or load a specific vignette
vignette("XQTL2_workflow", package = "XQTL2.Xplore")
vignette("XQTL2_usage", package = "XQTL2.Xplore")
```

## Troubleshooting

### Common Issues

1. **Package dependencies not found**
   ```r
   # Install missing dependencies
   install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork"))
   ```

2. **Bioconductor packages not available**
   ```r
   # Install BiocManager and required packages
   if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   BiocManager::install("vcfR")
   ```

3. **GitHub installation fails**
   - Make sure you have `devtools` or `remotes` installed
   - Check your internet connection
   - Verify the repository URL is correct

4. **Vignettes not loading**
   ```r
   # Reinstall with vignettes (this is the key fix!)
   devtools::install_github("tdlong/XQTL2.Xplore", build_vignettes = TRUE)
   ```

### System Requirements

- R version 4.0.0 or higher
- RStudio (recommended) or R console
- Internet connection for GitHub installation

## Data Requirements

The package includes example datasets, but for your own data, ensure it follows these formats:

- **QTL scan data**: Columns `chr`, `pos`, `Wald_log10p`, and optionally `cM`
- **Frequency data**: Columns `chr`, `pos`, `TRT`, `REP`, `founder`, `freq`
- **Gene data**: Columns `chr`, `start`, `end`, `gene_name`, `strand`, `is_utr`
- **Variant data**: Columns `CHROM`, `POS`, `type`, `subtype`, and genotype columns for each founder

## Support

If you encounter issues:

1. Check the vignettes for usage examples
2. Review the function documentation: `?XQTL_Manhattan_5panel`
3. Check the GitHub repository for updates
4. Open an issue on GitHub with a reproducible example

## Citation

If you use this package in your research, please cite:

```
XQTL2.Xplore: An R package for XQTL analysis and visualization
Author: Tony Long
Version: 0.0.0.9000
URL: https://github.com/tdlong/XQTL2.Xplore
``` 