#!/bin/bash
#SBATCH --job-name=find_common_genes
#SBATCH --output=logs_common_genes/%x_%j.out
#SBATCH --error=logs_common_genes/%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=64
#SBATCH --constraint=intel
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=srsanjeev@wisc.edu

eval "$(spack env activate --sh python_env)"
# Set the Python binary from spack environment
PYTHON_BIN="/mnt/scratch/home/wisc/srsanjeev/spack/var/spack/environments/python_env/.spack-env/view/bin/python"


# Input arguments
H5AD_DIR="patient_data_v2"
GENE_CSV="data/anal_precancerN5_cancerN21_hvg5k.csv"
OUTPUT_TXT="output/common_gene_n26_hvg5K.txt"

# Make sure logs & output directories exist
mkdir -p logs_common_genes
mkdir -p output

# Run the Python script
$PYTHON_PATH ./find_common_genes.py "$H5AD_DIR" "$GENE_CSV" "$OUTPUT_TXT"
