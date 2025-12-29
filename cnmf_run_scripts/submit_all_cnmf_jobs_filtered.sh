#!/bin/bash

# Master script to submit cNMF jobs for all datasets in a folder
# This version filters each dataset to common genes before running cNMF
# Usage: bash submit_all_cnmf_jobs_filtered.sh

# Create logs directory if it doesn't exist
mkdir -p logs_run_parallel_filtered

# Set the Python binary from spack environment
PYTHON_BIN="/home/wisc/srsanjeev/.conda/envs/amd_env/bin/python"

# Configuration - EDIT THESE PATHS
DATASET_FOLDER="patient_data_v2"
OUTPUT_DIR="cnmf_run-parallel_filtered_results" 
COMMON_GENES_FILE="output/common_gene_n26_hvg5K.txt"  # Path to common genes file
SBATCH_SCRIPT="run_cnmf_single_filtered.sh"        # The sbatch script in the current directory

# Job configuration
TOTAL_WORKERS=512      # Number of parallel workers
N_ITER=200          # Number of iterations
SEED=42             # Random seed

echo "====================================="
echo "Submitting cNMF jobs for all datasets"
echo "Dataset folder: $DATASET_FOLDER"
echo "Output directory: $OUTPUT_DIR"
echo "Common genes file: $COMMON_GENES_FILE"
echo "====================================="

# Check if common genes file exists
if [ ! -f "$COMMON_GENES_FILE" ]; then
    echo "ERROR: Common genes file not found at $COMMON_GENES_FILE"
    echo "Please ensure you have generated the common_genes.txt file"
    exit 1
fi

# Counter for submitted jobs
job_count=0

# Loop through all .h5ad files in the dataset folder
for counts_file in ${DATASET_FOLDER}/*.h5ad; do
    # Check if file exists (in case no files match the pattern)
    if [ ! -f "$counts_file" ]; then
        continue
    fi
    
    # Extract dataset name (without path and extension)
    dataset_name=$(basename "$counts_file" | sed 's/\.[^.]*$//')
    
    echo "Submitting job for dataset: $dataset_name"
    
    # Submit the job
    sbatch --job-name="cnmf_filtered_${dataset_name}" \
           --output="logs_run_parallel_filtered/cnmf_${dataset_name}_%j.out" \
           --error="logs_run_parallel_filtered/cnmf_${dataset_name}_%j.err" \
           ${SBATCH_SCRIPT} "$counts_file"
    
    ((job_count++))
done

echo "====================================="
echo "Submitted $job_count jobs"
echo "====================================="

