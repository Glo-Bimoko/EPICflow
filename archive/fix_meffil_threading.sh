#!/bin/bash
# fix_meffil_threading.sh
# Fix preprocessCore threading issues for meffil
# Based on: https://github.com/perishky/meffil/wiki/Common-problems

echo "============================================================"
echo "FIXING preprocessCore THREADING ISSUES FOR MEFFIL"
echo "============================================================"
echo ""
echo "This script implements the fix from meffil documentation:"
echo "https://github.com/perishky/meffil/wiki/Common-problems"
echo ""

# Activate conda environment
echo "Activating meffil_epicv2 environment..."
source ~/miniforge3/etc/profile.d/conda.sh
conda activate meffil_epicv2

echo ""
echo "Step 1: Removing existing preprocessCore..."
R --vanilla --quiet --no-save <<EOF
if ("preprocessCore" %in% rownames(installed.packages())) {
    remove.packages("preprocessCore")
    cat("✓ Removed old preprocessCore\n")
} else {
    cat("preprocessCore not currently installed\n")
}
EOF

echo ""
echo "Step 2: Installing preprocessCore WITHOUT multi-threading..."
echo "This is the recommended fix from meffil documentation"
echo ""

# The key fix: install with --disable-threading flag
R --vanilla --quiet --no-save <<EOF
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

cat("Installing preprocessCore with threading disabled...\n")
BiocManager::install("preprocessCore", 
                     configure.args = "--disable-threading",
                     update = FALSE,
                     ask = FALSE,
                     force = TRUE)

if (requireNamespace("preprocessCore", quietly = TRUE)) {
    cat("\n✓ preprocessCore installed successfully!\n")
    cat("Version:", as.character(packageVersion("preprocessCore")), "\n")
} else {
    cat("\n✗ Installation failed!\n")
    quit(status = 1)
}
EOF

if [ $? -ne 0 ]; then
    echo ""
    echo "Installation failed. Trying alternative method..."
    echo ""
    
    # Alternative: try conda installation
    echo "Attempting conda installation..."
    conda remove r-preprocesscore --force -y 2>/dev/null
    conda install -c bioconda r-preprocesscore -y
fi

echo ""
echo "Step 3: Testing the fix..."
R --vanilla --quiet --no-save <<EOF
library(preprocessCore)

# Test the problematic function
test_matrix <- matrix(rnorm(1000), ncol=10)
target <- rowMeans(test_matrix)

# Set single-threaded mode
Sys.setenv(OMP_NUM_THREADS = 1)

tryCatch({
    result <- normalize.quantiles.use.target(test_matrix, target)
    cat("\n✓ SUCCESS! preprocessCore is working correctly.\n")
    cat("The pthread error should be resolved.\n")
}, error = function(e) {
    cat("\n✗ Test failed:", conditionMessage(e), "\n")
    cat("Additional troubleshooting may be needed.\n")
    quit(status = 1)
})
EOF

if [ $? -eq 0 ]; then
    echo ""
    echo "============================================================"
    echo "FIX COMPLETE!"
    echo "============================================================"
    echo ""
    echo "You can now run your pipeline with functional normalization:"
    echo "  nextflow run main_meffil.nf -profile local \\"
    echo "    --idat_dir ./Plates_2_and_4/ \\"
    echo "    --out_dir ./results"
    echo ""
else
    echo ""
    echo "============================================================"
    echo "FIX INCOMPLETE"
    echo "============================================================"
    echo ""
    echo "Manual steps required:"
    echo "1. Start R in your conda environment"
    echo "2. Run these commands:"
    echo "   remove.packages('preprocessCore')"
    echo "   BiocManager::install('preprocessCore', configure.args='--disable-threading')"
    echo ""
fi
