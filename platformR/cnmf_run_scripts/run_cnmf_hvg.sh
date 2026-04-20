#!/bin/bash
#SBATCH --job-name=cnmf_hvg
#SBATCH --output=logs_cnmf_hvg/cnmf_%j.out
#SBATCH --error=logs_cnmf_hvg/cnmf_%j.err
#SBATCH --time=128:00:00
#SBATCH --mem=950G
#SBATCH --cpus-per-task=512
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=srsanjeev@wisc.edu

mkdir -p logs_cnmf_hvg

eval "$(conda shell.bash hook)"
conda activate amd_env
PYTHON_BIN="/home/wisc/srsanjeev/.conda/envs/amd_env/bin/python"

# Check if GNU parallel is available
if ! command -v parallel &> /dev/null; then
    echo "WARNING: GNU parallel not found. You may need to:"
    echo "  1. Load the module: module load parallel"
    echo "  2. Or install it: conda install -c conda-forge parallel"
fi

# ===== CONFIGURATION - EDIT THESE =====
COUNTS_FILE="path/to/your_combined.h5ad"   # <-- SET THIS to your .h5ad file
DATASET_NAME="combined_26patient"           # <-- SET THIS to a name for the run
OUTPUT_DIR="cnmf_hvg_results"
TOTAL_WORKERS=512
N_ITER=200
NUM_GENES=5000
SEED=14
# =======================================

echo "====================================="
echo "Running cNMF with HVG selection"
echo "Dataset: $DATASET_NAME"
echo "Counts file: $COUNTS_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "Num HVGs: $NUM_GENES"
echo "Workers: $TOTAL_WORKERS"
echo "Iterations: $N_ITER"
echo "====================================="

${PYTHON_BIN} ./run_parallel_hvg.py \
    --name ${DATASET_NAME} \
    --output-dir ${OUTPUT_DIR} \
    --counts ${COUNTS_FILE} \
    --numgenes ${NUM_GENES} \
    -k 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 \
    --n-iter ${N_ITER} \
    --total-workers ${TOTAL_WORKERS} \
    --seed ${SEED}

echo "Completed: $DATASET_NAME"
