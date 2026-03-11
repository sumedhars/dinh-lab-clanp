#!/bin/bash
#SBATCH --job-name=cnmf_filtered_single
#SBATCH --output=logs_run_parallel_filtered/cnmf_batch_%A_%a.out
#SBATCH --error=logs_run_parallel_filtered/cnmf_batch_%A_%a.err
#SBATCH --time=128:00:00
#SBATCH --mem=950G
#SBATCH --cpus-per-task=512
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=srsanjeev@wisc.edu

mkdir -p logs_run_parallel_filtered

# eval "$(spack env activate --sh python_env)"
eval "$(conda shell.bash hook)"
conda activate amd_env
PYTHON_BIN="/home/wisc/srsanjeev/.conda/envs/amd_env/bin/python"

# Check if GNU parallel is available
if ! command -v parallel &> /dev/null; then
    echo "WARNING: GNU parallel not found. You may need to:"
    echo "  1. Load the module: module load parallel"
    echo "  2. Or install it: conda install -c conda-forge parallel"
    echo "  3. Or use system package manager"
fi

# Configuration
DATASET_FOLDER="patient_data_v2"
OUTPUT_DIR="cnmf_run-parallel_filtered_results" 
COMMON_GENES_FILE="output/common_gene_n26_hvg5K.txt" 
TOTAL_WORKERS=512
N_ITER=200
SEED=14

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
echo "Common genes file: $COMMON_GENES_FILE"
echo "Python binary: $PYTHON_BIN"
echo "====================================="

# Check if common genes file exists
if [ ! -f "$COMMON_GENES_FILE" ]; then
    echo "ERROR: Common genes file not found at $COMMON_GENES_FILE"
    exit 1
fi

# Run cNMF with filtering to common genes
# Note: --numgenes is automatically set to the number of common genes
${PYTHON_BIN} ./run_parallel_filtered.py \
    --name ${DATASET_NAME} \
    --output-dir ${OUTPUT_DIR} \
    --counts ${COUNTS_FILE} \
    --common-genes-file ${COMMON_GENES_FILE} \
    -k 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 \
    --n-iter ${N_ITER} \
    --total-workers ${TOTAL_WORKERS} \
    --seed ${SEED}

echo "Completed processing dataset: $DATASET_NAME"
