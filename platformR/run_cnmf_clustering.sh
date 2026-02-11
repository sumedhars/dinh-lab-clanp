#!/bin/bash
#SBATCH --job-name=cnmf_clustering
#SBATCH --output=logs_consensus_clustering/cnmf_batch_%A_%a.out
#SBATCH --error=logs_consensus_clustering/cnmf_batch_%A_%a.err
#SBATCH --time=36:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=64
#SBATCH --constraint=intel
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=srsanjeev@wisc.edu

# ==============================================================================
# cNMF CONSENSUS + HIERARCHICAL CLUSTERING ACROSS DATASETS
# ==============================================================================
# This script performs the following workflow:
# 1. Runs cNMF consensus for each dataset with specified k values
# 2. Loads the consensus matrices (GEP x gene) from each dataset
# 3. Stacks all GEP matrices vertically into a unified matrix
# 4. Performs agglomerative hierarchical clustering on the GEPs
# 5. Creates visualizations to identify shared GEPs across patients/datasets
#
# Each dataset represents one patient, so this analysis reveals:
# - Whether patients share common gene expression programs (GEPs)
# - Which patients have similar transcriptional profiles
# - Patient-specific vs shared biological programs
# ==============================================================================

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================

# Path to Python binary in Spack environment
# This ensures we use the correct Python with all required packages
PYTHON_BIN="/mnt/scratch/home/wisc/srsanjeev/spack/var/spack/environments/python_env/.spack-env/view/bin/python"

# Create logs directory if it doesn't exist
# SLURM will write job output and errors here
mkdir -p logs_consensus_clustering

# Print job information header
echo "====================================================="
echo "cNMF Consensus + Hierarchical Clustering Analysis"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "Node: $SLURM_NODELIST"
echo "====================================================="

# ==============================================================================
# CONFIGURATION - ONLY TWO THINGS TO SET
# ==============================================================================

# 1. Set your parent folder path (contains dataset1, dataset2, etc.)
#    Each subdirectory should be a complete cNMF run directory with:
#    - prepared data from 'cnmf prepare'
#    - factorization results from 'cnmf factorize'
PARENT_DIR="cnmf_alg_results_subset"

# 2. Set K values for each dataset (in alphabetical order of dataset names)
#    These should be the optimal k values you determined from k_selection plots
#    Order: alphabetical by dataset folder name
SELECTED_KS="#, #, #, #, #, #, 11, 13, #, #, #, #, #, 15, #, #, #, #, #, #, #, #, 13, 8, #, #"

# ==============================================================================
# ADVANCED PARAMETERS (Optional - defaults usually work well)
# ==============================================================================

# Density threshold for cNMF consensus
# Lower values = more stringent, fewer but higher-confidence GEPs
# Higher values = more permissive, more GEPs included
# Default: 0.01 is standard for most analyses
DENSITY_THRESHOLD=0.01

# Hierarchical clustering method
# Options: 'average', 'complete', 'single', 'ward', 'weighted'
# - 'average' (UPGMA): Good balance, recommended for expression data
# - 'complete': Creates tight, compact clusters
# - 'ward': Minimizes within-cluster variance, good for clear groups
# - 'single': Can create chain-like clusters, less commonly used
LINKAGE_METHOD="average"

# Distance metric for clustering
# Options: 'correlation', 'euclidean', 'cosine', 'cityblock'
# - 'correlation': 1 - Pearson correlation, standard for expression data
# - 'euclidean': Standard distance in gene expression space
# - 'cosine': Focuses on pattern, ignores magnitude
# - 'cityblock': Manhattan distance, robust to outliers
DISTANCE_METRIC="correlation"

# Output prefix for all result files
# All output files will start with this prefix
OUTPUT_PREFIX="patient_gep_clustering"

# ==============================================================================
# PERFORMANCE OPTIMIZATION
# ==============================================================================

# Set number of threads for numerical libraries
# This optimizes performance for matrix operations
# Using all available CPUs allocated by SLURM
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo ""
echo "Configuration:"
echo "  Parent directory: ${PARENT_DIR}"
echo "  K values: ${SELECTED_KS}"
echo "  Density threshold: ${DENSITY_THRESHOLD}"
echo "  Linkage method: ${LINKAGE_METHOD}"
echo "  Distance metric: ${DISTANCE_METRIC}"
echo "  CPU threads: ${SLURM_CPUS_PER_TASK}"
echo "====================================================="

# ==============================================================================
# RUN ANALYSIS
# ==============================================================================

echo ""
echo "Starting analysis..."
echo ""

# Execute the Python script with all parameters
${PYTHON_BIN} cnmf_consensus_hierarchical_clustering.py \
    --parent_dir "${PARENT_DIR}" \
    --selected_ks "${SELECTED_KS}" \
    --density_threshold ${DENSITY_THRESHOLD} \
    --linkage_method "${LINKAGE_METHOD}" \
    --distance_metric "${DISTANCE_METRIC}" \
    --output_prefix "${OUTPUT_PREFIX}"

# ==============================================================================
# CHECK RESULTS AND REPORT
# ==============================================================================

# Check if the analysis completed successfully
if [ $? -eq 0 ]; then
    echo ""
    echo "====================================================="
    echo "ANALYSIS COMPLETED SUCCESSFULLY!"
    echo "====================================================="
    echo ""
    echo "Output files generated:"
    echo ""
    echo "1. ${OUTPUT_PREFIX}_stacked_gep_matrix.csv"
    echo "   → Full matrix of all GEPs from all patients"
    echo "   → Rows = GEPs (labeled by patient), Columns = genes"
    echo "   → Can be used for downstream analyses"
    echo ""
    echo "2. ${OUTPUT_PREFIX}_dendrogram.png"
    echo "   → Hierarchical clustering tree visualization"
    echo "   → GEPs are colored by patient/dataset"
    echo "   → Height shows distance at which GEPs cluster together"
    echo "   → Look for: mixed-patient clusters = shared programs"
    echo ""
    echo "3. ${OUTPUT_PREFIX}_distance_heatmap.png"
    echo "   → Pairwise distance matrix between all GEPs"
    echo "   → Dark regions = similar GEPs"
    echo "   → Off-diagonal blocks = cross-patient similarities"
    echo "   → Ordered by hierarchical clustering"
    echo ""
    echo "4. ${OUTPUT_PREFIX}_cluster_assignments.csv"
    echo "   → Table assigning each GEP to a cluster"
    echo "   → Columns: GEP name, Dataset, Cluster ID"
    echo "   → Use to identify which GEPs group together"
    echo ""
    echo "5. ${OUTPUT_PREFIX}_analysis_summary.txt"
    echo "   → Comprehensive text summary of results"
    echo "   → Dataset statistics, cluster composition"
    echo "   → Identifies shared vs patient-specific clusters"
    echo "   → Human-readable interpretation guide"
    echo ""
    echo "====================================================="
    echo "How to interpret results:"
    echo "====================================================="
    echo ""
    echo "SHARED PROGRAMS (cross-patient similarity):"
    echo "  • Dendrogram: Look for branches mixing different colored labels"
    echo "  • Heatmap: Look for dark off-diagonal blocks"
    echo "  • Clusters: Clusters with GEPs from multiple patients"
    echo "  → Suggests common biological programs across patients"
    echo ""
    echo "PATIENT-SPECIFIC PROGRAMS:"
    echo "  • Dendrogram: Branches with only one color/patient"
    echo "  • Heatmap: Dark blocks along the diagonal"
    echo "  • Clusters: Clusters dominated by one patient"
    echo "  → Suggests unique biology in that patient"
    echo ""
    echo "PATIENT SIMILARITY:"
    echo "  • Which patients' GEPs cluster together most often?"
    echo "  • Check summary file for cluster composition statistics"
    echo "  → Patients with similar transcriptional programs"
    echo ""
    echo "====================================================="
    echo "End time: $(date)"
    echo "====================================================="
else
    echo ""
    echo "====================================================="
    echo "ERROR: Analysis failed!"
    echo "====================================================="
    echo ""
    echo "Troubleshooting steps:"
    echo ""
    echo "1. Check the error log:"
    echo "   less logs/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err"
    echo ""
    echo "2. Verify your configuration:"
    echo "   • Does PARENT_DIR exist and contain dataset subdirectories?"
    echo "   • Do the subdirectories have cNMF prepare/factorize results?"
    echo "   • Does the number of K values match the number of datasets?"
    echo ""
    echo "3. Check individual dataset directories:"
    echo "   • Each should have been prepared with 'cnmf prepare'"
    echo "   • Each should have factorization results from 'cnmf factorize'"
    echo "   • Check for .cnmf_params.df.npz and factorize outputs"
    echo ""
    echo "4. Verify Python environment:"
    echo "   • Is the PYTHON_BIN path correct?"
    echo "   • Are required packages installed? (scipy, pandas, seaborn, cnmf)"
    echo ""
    echo "5. Check log files for detailed error messages"
    echo ""
    echo "====================================================="
    exit 1
fi