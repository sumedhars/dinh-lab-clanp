#!/bin/bash

# Master script to submit cNMF jobs for all datasets in a folder
# Usage: bash submit_all_cnmf_jobs.sh

# Set the Python binary from spack environment
PYTHON_BIN="/mnt/scratch/home/wisc/srsanjeev/spack/var/spack/environments/python_env/.spack-env/view/bin/python"

# Configuration 
DATASET_FOLDER="patient_data"
OUTPUT_DIR="cnmf_run-parallel_results" # Output directory for cNMF results
SBATCH_SCRIPT="run_cnmf_single.sh"

# Job configuration
TOTAL_WORKERS=64      # Number of parallel workers
N_ITER=200          # Number of iterations
SEED=14             # Random seed


echo "====================================="
echo "Submitting cNMF jobs for all datasets"
echo "Dataset folder: $DATASET_FOLDER"
echo "Output directory: $OUTPUT_DIR"
echo "====================================="

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
    sbatch --job-name="cnmf_std_${dataset_name}" \
           --output="logs_run_parallel/cnmf_${dataset_name}_%j.out" \
           --error="logs_run_parallel/cnmf_${dataset_name}_%j.err" \
           ${SBATCH_SCRIPT} "$counts_file"
    
    ((job_count++))
done

echo "====================================="
echo "Submitted $job_count jobs"
echo "====================================="
