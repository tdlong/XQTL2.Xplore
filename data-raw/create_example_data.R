#' Create example data objects from ZINC_Hanson files
#' 
#' This script loads the ZINC_Hanson example datasets and saves them as R data objects
#' for inclusion in the package. These datasets provide working examples for users.
#' 
#' Run this script when you need to update the example data.

library(usethis)
library(tibble)

# Read the pseudoscan data (QTL scan results)
zinc_hanson_pseudoscan <- as_tibble(read.table("../bitnbobs/ZINC_Hanson.pseudoscan.txt"))

# Read the means by sample data (frequency data)
zinc_hanson_means <- as_tibble(read.table("../bitnbobs/ZINC_Hanson.meansBySample.txt"))

# Save the data objects
usethis::use_data(zinc_hanson_pseudoscan, overwrite = TRUE)
usethis::use_data(zinc_hanson_means, overwrite = TRUE)

# Print summary information
cat("ZINC Hanson pseudoscan data:\n")
cat("  Dimensions:", dim(zinc_hanson_pseudoscan), "\n")
cat("  Columns:", paste(names(zinc_hanson_pseudoscan), collapse = ", "), "\n")
cat("  Chromosomes:", paste(unique(zinc_hanson_pseudoscan$chr), collapse = ", "), "\n")
cat("  Position range:", range(zinc_hanson_pseudoscan$pos, na.rm = TRUE), "\n\n")

cat("ZINC Hanson means data:\n")
cat("  Dimensions:", dim(zinc_hanson_means), "\n")
cat("  Columns:", paste(names(zinc_hanson_means), collapse = ", "), "\n")
cat("  Chromosomes:", paste(unique(zinc_hanson_means$chr), collapse = ", "), "\n")
cat("  Position range:", range(zinc_hanson_means$pos, na.rm = TRUE), "\n")
cat("  Treatments:", paste(unique(zinc_hanson_means$TRT), collapse = ", "), "\n")
cat("  Replicates:", paste(unique(zinc_hanson_means$REP), collapse = ", "), "\n")
cat("  Founders:", paste(unique(zinc_hanson_means$founder), collapse = ", "), "\n") 