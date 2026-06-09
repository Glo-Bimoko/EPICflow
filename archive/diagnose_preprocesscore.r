#!/usr/bin/env Rscript
# fix_preprocesscore.R
# Reinstall preprocessCore without threading support

cat("\n", rep("=", 60), "\n", sep = "")
cat("FIXING preprocessCore THREADING ISSUES\n")
cat(rep("=", 60), "\n\n")

cat("This script will reinstall preprocessCore without OpenMP threading\n")
cat("to fix the pthread_create error 22\n\n")

# Check if BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

library(BiocManager)

# Remove existing preprocessCore
cat("Step 1: Removing existing preprocessCore...\n")
tryCatch({
  remove.packages("preprocessCore")
  cat("✓ Removed old preprocessCore\n\n")
}, error = function(e) {
  cat("Note: preprocessCore was not installed or couldn't be removed\n\n")
})

# Install preprocessCore with threading disabled
cat("Step 2: Installing preprocessCore WITHOUT threading...\n")
cat("This may take a few minutes...\n\n")

# Method 1: Try with configure.args (best option)
success <- FALSE
tryCatch({
  BiocManager::install(
    "preprocessCore",
    configure.args = "--disable-threading",
    update = FALSE,
    ask = FALSE,
    force = TRUE
  )
  success <- TRUE
  cat("✓ Successfully installed preprocessCore without threading!\n\n")
}, error = function(e) {
  cat("✗ Method 1 failed:", conditionMessage(e), "\n")
  cat("Trying alternative method...\n\n")
})

# Method 2: If Method 1 fails, try forcing source installation
if (!success) {
  tryCatch({
    # Set environment to disable OpenMP before installation
    Sys.setenv(MAKEFLAGS = "-j1")
    Sys.setenv(R_COMPILE_PKGS = "1")
    
    BiocManager::install(
      "preprocessCore",
      type = "source",
      update = FALSE,
      ask = FALSE,
      force = TRUE
    )
    success <- TRUE
    cat("✓ Successfully installed preprocessCore from source!\n\n")
  }, error = function(e) {
    cat("✗ Method 2 also failed:", conditionMessage(e), "\n\n")
  })
}

# Verify installation
cat("Step 3: Verifying installation...\n")
if (requireNamespace("preprocessCore", quietly = TRUE)) {
  library(preprocessCore)
  cat("✓ preprocessCore version:", packageVersion("preprocessCore"), "\n")
  
  # Quick test
  cat("\nStep 4: Testing with sample data...\n")
  test_matrix <- matrix(rnorm(100), ncol=10)
  target <- rowMeans(test_matrix)
  
  tryCatch({
    Sys.setenv(OMP_NUM_THREADS = 1)
    result <- preprocessCore::normalize.quantiles.use.target(test_matrix, target)
    cat("✓ Test passed! preprocessCore is working correctly!\n\n")
    
    cat(rep("=", 60), "\n")
    cat("SUCCESS! You can now run your pipeline.\n")
    cat(rep("=", 60), "\n\n")
    
  }, error = function(e) {
    cat("✗ Test failed:", conditionMessage(e), "\n")
    cat("\nThe installation succeeded but threading issues persist.\n")
    cat("You may need to use the alternative normalization approach.\n\n")
  })
  
} else {
  cat("✗ preprocessCore installation failed!\n")
  cat("\nTrying one more approach: conda installation...\n")
  cat("Run this in your terminal:\n")
  cat("  conda activate meffil_epicv2\n")
  cat("  conda remove r-preprocesscore\n")
  cat("  conda install -c conda-forge r-preprocesscore\n\n")
}
