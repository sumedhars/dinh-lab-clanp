#!/bin/bash
#SBATCH --job-name=cnmf_batch
#SBATCH --output=logs_run_parallel/cnmf_batch_%A_%a.out
#SBATCH --error=logs_run_parallel/cnmf_batch_%A_%a.err
#SBATCH --time=54:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=64
#SBATCH --constraint=intel
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=srsanjeev@wisc.edu

mkdir -p logs_run_parallel

eval "$(spack env activate --sh python_env)"
# Set the Python binary from spack environment
PYTHON_BIN="/mnt/scratch/home/wisc/srsanjeev/spack/var/spack/environments/python_env/.spack-env/view/bin/python"

# Check if GNU parallel is available
if ! command -v parallel &> /dev/null; then
    echo "WARNING: GNU parallel not found. You may need to:"
    echo "  1. Load the module: module load parallel"
    echo "  2. Or install it: conda install -c conda-forge parallel"
    echo "  3. Or use system package manager"
fi

# Configuration
DATASET_FOLDER="patient_data" 
OUTPUT_DIR="cnmf_run-parallel_results"  
TOTAL_WORKERS=64      # Number of parallel workers
N_ITER=200          # Number of iterations
SEED=14             # Random seed

# Get the dataset file from array task ID if using job arrays
# Otherwise, you can pass it as an argument
if [ ! -z "$SLURM_ARRAY_TASK_ID" ]; then
    # If running as job array, get dataset based on array index
    DATASETS=(${DATASET_FOLDER}/*.h5ad)
    COUNTS_FILE="${DATASETS[$SLURM_ARRAY_TASK_ID]}"
else
    # If running single job, use command line argument
    COUNTS_FILE=$1
fi

# Extract dataset name (without path and extension)
DATASET_NAME=$(basename "$COUNTS_FILE" | sed 's/\.[^.]*$//')

echo "====================================="
echo "Processing dataset: $DATASET_NAME"
echo "Counts file: $COUNTS_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "Python binary: $PYTHON_BIN"
echo "====================================="

# Run cNMF with k=5-21 and numgenes=2000
${PYTHON_BIN} ./run_parallel_modified.py \
    --name ${DATASET_NAME} \
    --output-dir ${OUTPUT_DIR} \
    --counts ${COUNTS_FILE} \
    -k 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 \
    --n-iter ${N_ITER} \
    --total-workers ${TOTAL_WORKERS} \
    --numgenes 2000 \
    --seed ${SEED}

echo "Completed processing dataset: $DATASET_NAME"


