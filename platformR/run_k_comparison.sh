#!/bin/bash
#SBATCH --job-name=k_comparison
#SBATCH --output=logs_k_comparison/%x_%j.out
#SBATCH --error=logs_k_comparison/%x_%j.err
#SBATCH --time=54:00:00
#SBATCH --mem=950G
#SBATCH --cpus-per-task=512
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=srsanjeev@wisc.edu

# ==============================================================================
# cNMF K-VALUE COMPARISON: HIERARCHICAL CLUSTERING ANALYSIS
# ==============================================================================
# This script analyzes a SINGLE dataset across MULTIPLE k values to identify:
#   - Which GEPs are consistently found across different k values
#   - How GEPs cluster together (hierarchical clustering)
#   - Correlation structure between GEPs from different k values
#
# This helps you determine:
#   - Whether certain biological programs are robust across k choices
#   - What k value best captures your biological processes of interest
#   - Which GEPs might be redundant (highly correlated within same k)
# ==============================================================================

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================

eval "$(conda shell.bash hook)"
conda activate amd_env
PYTHON_BIN="/home/wisc/srsanjeev/.conda/envs/amd_env/bin/python"

# Create logs directory if it doesn't exist
mkdir -p logs_k_comparison

# Print job information
echo "====================================================="
echo "cNMF K-Value Comparison Analysis"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "Node: $SLURM_NODELIST"
echo "====================================================="

# ==============================================================================
# CONFIGURATION - UPDATE THESE TWO SETTINGS
# ==============================================================================

# 1. Path to your cNMF results directory for a SINGLE dataset
CNMF_DIR="cnmf_run-parallel_filtered_results/anal_pc5_c21_S13.filtered"  # ← UPDATE THIS!

# 2. K values to compare (space-separated)
K_VALUES="5 6 7 8"  # ← UPDATE THIS!

# ==============================================================================
# ADVANCED PARAMETERS (Optional)
# ==============================================================================

# Density threshold for cNMF consensus
# This should match what you used in your main analysis
# Default: 0.01
DENSITY_THRESHOLD=0.01

# Hierarchical clustering linkage method
# Options: 'average', 'complete', 'single', 'ward', 'weighted'
# - 'average' (UPGMA): Good balance, recommended for expression data
# - 'ward': Minimizes within-cluster variance
# - 'complete': Creates tight, compact clusters
LINKAGE_METHOD="average"

# Output prefix for all result files
OUTPUT_PREFIX="S13_k5_k6_k7_k8_comparison"

# Skip consensus step if you already ran it
# Set to "--skip_consensus" if consensus files already exist
# Leave empty ("") to run consensus
SKIP_CONSENSUS=""

# ==============================================================================
# PERFORMANCE OPTIMIZATION
# ==============================================================================

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo ""
echo "Configuration:"
echo "  Dataset directory: ${CNMF_DIR}"
echo "  K values: ${K_VALUES}"
echo "  Density threshold: ${DENSITY_THRESHOLD}"
echo "  Linkage method: ${LINKAGE_METHOD}"
echo "  Output prefix: ${OUTPUT_PREFIX}"
echo "  CPU threads: ${SLURM_CPUS_PER_TASK}"
if [ -n "$SKIP_CONSENSUS" ]; then
    echo "  Skip consensus: YES (using existing files)"
else
    echo "  Skip consensus: NO (will run consensus)"
fi
echo "====================================================="

# ==============================================================================
# RUN ANALYSIS
# ==============================================================================

echo ""
echo "Starting analysis..."
echo ""

# Build command
CMD="${PYTHON_BIN} cnmf_k_comparison_clustering.py \
    --cnmf_dir ${CNMF_DIR} \
    --k_values ${K_VALUES} \
    --density_threshold ${DENSITY_THRESHOLD} \
    --linkage_method ${LINKAGE_METHOD} \
    --output_prefix ${OUTPUT_PREFIX}"

# Add skip consensus flag if set
if [ -n "$SKIP_CONSENSUS" ]; then
    CMD="${CMD} ${SKIP_CONSENSUS}"
fi

# Execute
eval $CMD

# ==============================================================================
# CHECK RESULTS AND REPORT
# ==============================================================================


if [ $? -eq 0 ]; then
    echo ""
    echo "====================================================="
    echo "ANALYSIS COMPLETED SUCCESSFULLY!"
    echo "====================================================="
    echo ""
    echo "Output files generated:"
    echo ""
    echo "1. ${OUTPUT_PREFIX}_stacked_gep_matrix.csv"
    echo "   → Full matrix of all GEPs from all k values"
    echo "   → Rows = GEPs (labeled as k8_GEP1, k10_GEP2, etc.)"
    echo "   → Columns = genes with Z-score weights"
    echo "   → Use for downstream analysis or manual inspection"
    echo ""
    echo "2. ${OUTPUT_PREFIX}_correlation_matrix.csv"
    echo "   → Pearson correlation between all pairs of GEPs"
    echo "   → Values range from -1 (anti-correlated) to +1 (perfectly correlated)"
    echo "   → High cross-k correlations suggest stable programs"
    echo ""
    echo "3. ${OUTPUT_PREFIX}_dendrogram.png"
    echo "   → Hierarchical clustering tree visualization"
    echo "   → GEPs are colored by k value"
    echo "   → Height shows distance (1 - correlation)"
    echo "   → Look for: mixed-color clusters = stable across k values"
    echo ""
    echo "4. ${OUTPUT_PREFIX}_correlation_heatmap.png"
    echo "   → Clustered heatmap of all pairwise correlations"
    echo "   → Red = positive correlation, Blue = negative correlation"
    echo "   → Color bars on axes show which k value each GEP comes from"
    echo "   → Diagonal blocks = GEPs from same k value"
    echo "   → Off-diagonal blocks = cross-k correlations"
    echo ""
    
    # Count number of k values to determine if pairwise files exist
    K_COUNT=$(echo ${K_VALUES} | wc -w)
    
    if [ ${K_COUNT} -le 4 ]; then
        echo "5. ${OUTPUT_PREFIX}_pairwise_k*_vs_k*.png (multiple files)"
        echo "   → Individual pairwise comparison heatmaps for each k-value pair"
        echo "   → Example: k4 vs k5, k4 vs k6, k5 vs k6 for k_values=\"4 5 6\""
        echo "   → Each heatmap shows correlations between GEPs from two k values"
        echo "   → Rows = GEPs from lower k, Columns = GEPs from higher k"
        echo "   → Values shown in cells for small matrices (≤12 GEPs per k)"
        echo "   → Generated only when comparing 4 or fewer k values"
        echo ""
        echo "6. ${OUTPUT_PREFIX}_similar_geps_threshold*.csv"
        echo "   → Tables of highly correlated GEP pairs"
        echo "   → Generated for thresholds: 0.9, 0.8, 0.7"
        echo "   → Includes both cross-k and same-k pairs"
        echo "   → Focus on cross-k pairs for stable programs"
        echo ""
        echo "7. ${OUTPUT_PREFIX}_analysis_summary.txt"
        echo "   → Comprehensive text summary with interpretation guide"
        echo "   → Statistics, correlation distributions"
        echo "   → Guidance on how to use results for k selection"
        echo "   → **READ THIS FIRST** for understanding your results"
    else
        echo "5. ${OUTPUT_PREFIX}_similar_geps_threshold*.csv"
        echo "   → Tables of highly correlated GEP pairs"
        echo "   → Generated for thresholds: 0.9, 0.8, 0.7"
        echo "   → Includes both cross-k and same-k pairs"
        echo "   → Focus on cross-k pairs for stable programs"
        echo ""
        echo "6. ${OUTPUT_PREFIX}_analysis_summary.txt"
        echo "   → Comprehensive text summary with interpretation guide"
        echo "   → Statistics, correlation distributions"
        echo "   → Guidance on how to use results for k selection"
        echo "   → **READ THIS FIRST** for understanding your results"
        echo ""
        echo "NOTE: Pairwise k-value comparison heatmaps are only generated"
        echo "      when comparing 4 or fewer k values. You are comparing ${K_COUNT} k values."
    fi
    echo ""
    echo "====================================================="
    echo "How to interpret results:"
    echo "====================================================="
    echo ""
    echo "CONSISTENT GEPs (across k values):"
    echo "  • Dendrogram: Look for clusters mixing different k values"
    echo "  • Heatmap: Look for red (high correlation) off-diagonal blocks"
    echo "  • Similar GEPs table: High cross-k correlations (≥0.8)"
    echo "  → These are robust biological programs"
    echo "  → Less sensitive to k parameter choice"
    echo ""
    echo "K-SPECIFIC GEPs:"
    echo "  • Dendrogram: Branches containing only one k value"
    echo "  • Heatmap: Strong correlations within diagonal blocks only"
    echo "  • Similar GEPs table: Low or no cross-k correlations"
    echo "  → May represent finer biological distinctions"
    echo "  → Or over-fitting / splitting of programs"
    echo ""
    echo "CHOOSING OPTIMAL K:"
    echo "  • Look for k where:"
    echo "    1. Key stable programs are captured (seen in multiple k values)"
    echo "    2. GEPs are distinct (not too redundant within k)"
    echo "    3. Biological resolution matches your research question"
    echo "  • Compare gene loadings for highly correlated cross-k GEPs"
    echo "  • Consider using the lower k that captures stable programs"
    echo "    (simpler model, easier interpretation)"
    echo ""
    echo "NEXT STEPS:"
    echo "  1. Read the analysis_summary.txt file for detailed interpretation"
    echo "  2. Examine dendrogram and heatmap for stable GEP clusters"
    echo "  3. Check similar_geps tables for cross-k correlations"
    echo "  4. For highly correlated GEPs, compare their top genes"
    echo "  5. Run gene set enrichment on stable programs"
    echo "  6. Select k based on biological insight + statistical stability"
    echo ""
    echo "====================================================="
    echo "End time: $(date)"
    echo "====================================================="
else
    echo ""
    echo "====================================================="
    echo "ERROR: Analysis failed!"
    echo "====================================================="
    exit 1
fi
