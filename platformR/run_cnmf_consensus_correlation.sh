#!/bin/bash
#SBATCH --job-name=cnmf_simple
#SBATCH --output=logs_corr/%x_%j.out
#SBATCH --error=logs_corr/%x_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=32
#SBATCH --constraint=intel
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=srsanjeev@wisc.edu

# ==============================================================================
# SIMPLE cNMF CONSENSUS AND CORRELATION - MINIMAL CONFIGURATION
# ==============================================================================
# This is the simplest version - just set your parent folder and K values!


# Activate Spack Python environment
spack env activate python_env
# Path to your Spack Python environment
PYTHON_BIN="/mnt/scratch/home/wisc/srsanjeev/spack/var/spack/environments/python_env/.spack-env/view/bin/python"

mkdir -p logs_corr

echo "====================================================="
echo "cNMF Consensus and Correlation Analysis"
echo "Job ID: $SLURM_JOB_ID"
echo "====================================================="

# ==============================================================================
# ONLY TWO THINGS TO CONFIGURE:
# ==============================================================================

# 1. Set your parent folder path (contains dataset1, dataset2, etc.)
PARENT_DIR="cnmf_alg_results_subset"  

# 2. Set K values for each dataset (in alphabetical order of dataset names)
#    For 4 datasets: dataset1, dataset2, dataset3, dataset4
SELECTED_KS="10,10,10,10"  # ← UPDATE THESE based on k_selection plots!

# ==============================================================================
# RUN ANALYSIS (no need to modify below)
# ==============================================================================

# Set CPU threads for optimal performance
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Run the analysis
${PYTHON_BIN} cnmf_gep_gene_corr.py \
    --parent_dir "${PARENT_DIR}" \
    --selected_ks "${SELECTED_KS}" \
    --density_threshold 0.01 \
    --output_prefix "cnmf_correlation_results" \
    --prefer_spectra tpm \
    --plot_global_heatmap \
    --draw_dataset_blocks
    
# Check if successful
if [ $? -eq 0 ]; then
    echo "====================================================="
    echo "Analysis Complete! Output files:"
    echo "  - cnmf_correlation_results_correlation_matrix.csv"
    echo "  - cnmf_correlation_results_p_values.csv"  
    echo "  - cnmf_correlation_results_correlation_heatmap.png"
    echo "  - cnmf_correlation_results_dataset_summary.csv"
    echo "====================================================="
else
    echo "ERROR: Analysis failed. Check logs for details."
    exit 1
fi
