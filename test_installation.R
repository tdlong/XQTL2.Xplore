# Test script to verify XQTL2.Xplore installation
# Run this after installation to make sure everything works

cat("Testing XQTL2.Xplore installation...\n\n")

# Test 1: Load the package
cat("1. Loading package...\n")
library(XQTL2.Xplore)
cat("   ✓ Package loaded successfully\n\n")

# Test 2: Load example datasets
cat("2. Loading example datasets...\n")
data(zinc_hanson_pseudoscan)
data(zinc_hanson_means)
cat("   ✓ Example datasets loaded\n")
cat("   - zinc_hanson_pseudoscan:", dim(zinc_hanson_pseudoscan), "\n")
cat("   - zinc_hanson_means:", dim(zinc_hanson_means), "\n\n")

# Test 3: Test basic function
cat("3. Testing basic function...\n")
p <- XQTL_Manhattan_5panel(zinc_hanson_pseudoscan, cM = FALSE)
cat("   ✓ Manhattan plot created successfully\n\n")

# Test 4: Check vignettes
cat("4. Checking vignettes...\n")
vignettes <- vignette(package = "XQTL2.Xplore")
cat("   ✓ Vignettes available:", paste(vignettes$results[, "Item"], collapse = ", "), "\n\n")

cat("🎉 Installation test completed successfully!\n")
cat("Your XQTL2.Xplore package is ready to use.\n")
cat("Run 'browseVignettes(\"XQTL2.Xplore\")' to see tutorials.\n") 