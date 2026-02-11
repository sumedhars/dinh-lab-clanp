#!/bin/bash
#SBATCH --job-name=consensus_runner
#SBATCH --output=logs_consensus_runner/%x_%j.out
#SBATCH --error=logs_consensus_runner/%x_%j.err
#SBATCH --time=20:00:00
#SBATCH --mem=450G
#SBATCH --cpus-per-task=150
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=srsanjeev@wisc.edu

# ==============================================================================
# cNMF CONSENSUS RUNNER: Run consensus for multiple k values
# ==============================================================================
# This script runs cNMF consensus analysis for multiple k values and organizes
# the output files into a single folder for easy access.
#
# Prerequisites:
#   - You must have already run: cnmf prepare, factorize, combine, k_selection_plot
#   - All k values you want to process must have completed factorization
#
# What this script does:
#   1. Runs 'cnmf consensus' for each k value in your list
#   2. Creates a folder with organized consensus outputs
#   3. Copies these files for each k value:
#      - gene_spectra_score (Z-score GEP matrix)
#      - gene_spectra_tpm (TPM GEP matrix)  
#      - usages (cell usage matrix)
#      - clustering plots (if requested)
# ==============================================================================

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================

eval "$(conda shell.bash hook)"
conda activate amd_env
PYTHON_BIN="/home/wisc/srsanjeev/.conda/envs/amd_env/bin/python"

# Create logs directory if it doesn't exist
mkdir -p logs_consensus_runner

# Print job information
echo "====================================================="
echo "cNMF Consensus Runner"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "Node: $SLURM_NODELIST"
echo "====================================================="

# ==============================================================================
# CONFIGURATION - UPDATE THESE SETTINGS
# ==============================================================================

# 1. Full path to your cNMF dataset directory
#    This should be the directory that contains your consensus outputs
#    Structure: /parent_dir/dataset_name/
#               /parent_dir/dataset_name/cnmf_tmp/
#    Example: /path/to/results/patient1
DATASET_DIR="cnmf_run-parallel_filtered_results/anal_pc5_c21_S18.filtered"  # ← UPDATE THIS!

# 2. K values to run consensus on (space-separated)
K_VALUES="5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21"  # ← UPDATE THIS!

# ==============================================================================
# OPTIONAL PARAMETERS
# ==============================================================================

# Density threshold for consensus
# This filters out low-quality programs
# Default: 0.01 (recommended)
# Higher values = stricter filtering
DENSITY_THRESHOLD=0.01

# Name of output folder to create (optional)
# Leave empty to use default: <DATASET_NAME>_consensus_outputs
OUTPUT_FOLDER=""

# Generate clustering diagnostic plots?
# Set to "--show_clustering" to generate plots
# Leave empty ("") to skip plots
SHOW_CLUSTERING=""

# ==============================================================================
# PERFORMANCE OPTIMIZATION
# ==============================================================================

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo ""
echo "Configuration:"
echo "  Dataset directory: ${DATASET_DIR}"
if [ -n "$OUTPUT_FOLDER" ]; then
    echo "  Output folder: ${OUTPUT_FOLDER}"
else
    DATASET_NAME=$(basename "${DATASET_DIR}")
    echo "  Output folder: ${DATASET_NAME}_consensus_outputs (default)"
fi
echo "  K values: ${K_VALUES}"
echo "  Density threshold: ${DENSITY_THRESHOLD}"
if [ -n "$SHOW_CLUSTERING" ]; then
    echo "  Clustering plots: YES"
else
    echo "  Clustering plots: NO"
fi
echo "  CPU threads: ${SLURM_CPUS_PER_TASK}"
echo "====================================================="

# ==============================================================================
# RUN CONSENSUS ANALYSIS
# ==============================================================================

echo ""
echo "Starting consensus analysis..."
echo ""

# Build command
CMD="${PYTHON_BIN} cnmf_consensus_runner.py \
    --dataset_dir ${DATASET_DIR} \
    --k_values ${K_VALUES} \
    --density_threshold ${DENSITY_THRESHOLD}"

# Add optional parameters if set
if [ -n "$OUTPUT_FOLDER" ]; then
    CMD="${CMD} --output_folder ${OUTPUT_FOLDER}"
fi

if [ -n "$SHOW_CLUSTERING" ]; then
    CMD="${CMD} ${SHOW_CLUSTERING}"
fi

# Execute
eval $CMD

# ==============================================================================
# CHECK RESULTS AND REPORT
# ==============================================================================

if [ $? -eq 0 ]; then
    echo ""
    echo "====================================================="
    echo "✓ CONSENSUS ANALYSIS COMPLETED SUCCESSFULLY!"
    echo "====================================================="
    echo ""
    
    # Extract dataset name from path
    DATASET_NAME=$(basename "${DATASET_DIR}")
    
    # Determine actual output folder name
    if [ -n "$OUTPUT_FOLDER" ]; then
        ACTUAL_OUTPUT="${OUTPUT_FOLDER}"
    else
        ACTUAL_OUTPUT="${DATASET_NAME}_consensus_outputs"
    fi
    
    echo "All consensus files are organized in:"
    echo "  ./${ACTUAL_OUTPUT}/"
    echo ""
    echo "Files generated for each k value:"
    echo "  1. ${DATASET_NAME}.gene_spectra_score.k_<k>.dt_${DENSITY_THRESHOLD//./_}.txt"
    echo "     → Z-score normalized GEP matrix"
    echo "     → Use for identifying top genes in each program"
    echo ""
    echo "  2. ${DATASET_NAME}.gene_spectra_tpm.k_<k>.dt_${DENSITY_THRESHOLD//./_}.txt"
    echo "     → TPM-normalized GEP matrix"
    echo "     → Use for comparing expression levels"
    echo ""
    echo "  3. ${DATASET_NAME}.usages.k_<k>.dt_${DENSITY_THRESHOLD//./_}.consensus.txt"
    echo "     → Cell usage matrix (cells x GEPs)"
    echo "     → Each row sums to 1"
    echo "     → Use for downstream analysis of program activity"
    echo ""
    
    if [ -n "$SHOW_CLUSTERING" ]; then
        echo "  4. ${DATASET_NAME}.clustering.k_<k>.dt_${DENSITY_THRESHOLD//./_}.pdf"
        echo "     → Diagnostic clustergram showing GEP relationships"
        echo ""
    fi
    
    echo "A summary file has been created:"
    echo "  ./${ACTUAL_OUTPUT}/consensus_run_summary.txt"
    echo ""
    echo "====================================================="
    echo "NEXT STEPS"
    echo "====================================================="
    echo ""
    echo "1. Review the summary file for run details"
    echo ""
    echo "2. Load and analyze your GEP matrices:"
    echo "   import pandas as pd"
    echo "   gep_scores = pd.read_csv('./${ACTUAL_OUTPUT}/${DATASET_NAME}.gene_spectra_score.k_10.dt_0_01.txt',"
    echo "                            sep='\t', index_col=0)"
    echo ""
    echo "3. Load usage matrix to see GEP activity:"
    echo "   usage = pd.read_csv('./${ACTUAL_OUTPUT}/${DATASET_NAME}.usages.k_10.dt_0_01.consensus.txt',"
    echo "                       sep='\t', index_col=0)"
    echo ""
    echo "4. Compare GEPs across k values:"
    echo "   Use cnmf_k_comparison_clustering.py to identify stable programs"
    echo ""
    echo "5. Perform gene set enrichment:"
    echo "   - Extract top genes for each GEP (high Z-scores)"
    echo "   - Run GO/KEGG/Reactome enrichment"
    echo "   - Interpret biological meaning of programs"
    echo ""
    echo "6. Visualize GEP usage patterns:"
    echo "   - Plot usage across cell types"
    echo "   - Create UMAP colored by GEP usage"
    echo "   - Correlate with cell metadata"
    echo ""
    echo "====================================================="
    echo "End time: $(date)"
    echo "====================================================="
else
    echo ""
    echo "====================================================="
    echo "❌ ERROR: Consensus analysis failed!"
    echo "====================================================="
    echo ""
    exit 1
fi
