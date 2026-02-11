#!/usr/bin/env python3
"""
cNMF Consensus with Hierarchical Clustering Across Datasets
============================================================

This script performs the following workflow:
1. Runs cNMF consensus for each dataset with specified k values
2. Loads the consensus matrices (GEP x gene) from each dataset
3. Stacks all GEP matrices vertically into a unified matrix
4. Performs agglomerative hierarchical clustering on GEPs
5. Creates visualizations to identify shared GEPs across patients/datasets

Author: Sumedha
Date: November 2025
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr
import subprocess
import warnings
warnings.filterwarnings('ignore')

# Set plotting style
sns.set_style("white")
plt.rcParams['font.size'] = 10
plt.rcParams['figure.dpi'] = 300


def parse_arguments():
    """
    Parse command line arguments.
    
    Returns:
        argparse.Namespace: Parsed arguments containing:
            - parent_dir: Path to parent directory containing dataset subdirectories
            - selected_ks: Comma-separated k values for each dataset
            - density_threshold: Threshold for cNMF consensus (default: 0.01)
            - output_prefix: Prefix for output files
            - linkage_method: Method for hierarchical clustering
            - distance_metric: Distance metric for clustering
    """
    parser = argparse.ArgumentParser(
        description='Run cNMF consensus and hierarchical clustering across datasets',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python cnmf_consensus_hierarchical_clustering.py \\
        --parent_dir /path/to/parent_folder \\
        --selected_ks "10,12,11,10" \\
        --output_prefix patient_gep_clustering
        """
    )
    
    parser.add_argument('--parent_dir', type=str, required=True,
                        help='Parent directory containing dataset subdirectories')
    parser.add_argument('--selected_ks', type=str, required=True,
                        help='Comma-separated k values (one per dataset, in alphabetical order)')
    parser.add_argument('--density_threshold', type=float, default=0.01,
                        help='Density threshold for cNMF consensus (default: 0.01)')
    parser.add_argument('--output_prefix', type=str, default='cnmf_clustering',
                        help='Prefix for output files (default: cnmf_clustering)')
    parser.add_argument('--linkage_method', type=str, default='average',
                        choices=['single', 'complete', 'average', 'weighted', 'ward'],
                        help='Linkage method for hierarchical clustering (default: average)')
    parser.add_argument('--distance_metric', type=str, default='correlation',
                        choices=['correlation', 'euclidean', 'cosine', 'cityblock'],
                        help='Distance metric for clustering (default: correlation)')
    
    return parser.parse_args()


def format_density_threshold(threshold):
    """
    Format density threshold for cNMF filenames.
    cNMF replaces the decimal point with an underscore in output filenames.
    For example: 0.01 becomes 0_01, 0.1 becomes 0_1
    
    Args:
        threshold (float): Density threshold value
        
    Returns:
        str: Formatted threshold string for filenames
    """
    return str(threshold).replace('.', '_')


def get_dataset_directories(parent_dir):
    """
    Get sorted list of dataset subdirectories.
    
    Args:
        parent_dir (str): Path to parent directory
        
    Returns:
        list: Sorted list of dataset directory names
        
    Purpose: Identify all dataset subdirectories in alphabetical order
    to match with the k values provided in the same order.
    """
    parent_path = Path(parent_dir)
    if not parent_path.exists():
        raise ValueError(f"Parent directory does not exist: {parent_dir}")
    
    # Get all subdirectories
    datasets = sorted([d.name for d in parent_path.iterdir() if d.is_dir()])
    
    if len(datasets) == 0:
        raise ValueError(f"No subdirectories found in {parent_dir}")
    
    print(f"\nFound {len(datasets)} datasets:")
    for i, dataset in enumerate(datasets, 1):
        print(f"  {i}. {dataset}")
    
    return datasets


def run_cnmf_consensus(dataset_dir, k_value, density_threshold=0.01):
    """
    Run cNMF consensus command for a single dataset.
    
    Args:
        dataset_dir (Path): Path to dataset directory
        k_value (int): Number of gene expression programs (k)
        density_threshold (float): Density threshold for consensus
        
    Returns:
        Path: Path to the consensus spectra file (GEP x gene matrix)
        
    Purpose: Execute the cNMF consensus step which:
    - Combines multiple cNMF runs (typically from prepare/factorize steps)
    - Identifies consensus gene expression programs
    - Outputs a spectra file where rows=GEPs, columns=genes
    """
    print(f"\n{'='*70}")
    print(f"Running cNMF consensus for: {dataset_dir.name}")
    print(f"K value: {k_value}")
    print(f"{'='*70}")
    
    # Construct the consensus command
    # This follows the standard cNMF workflow: prepare -> factorize -> consensus -> usage
    # IMPORTANT: --output-dir must be the PARENT directory, --name is the subdirectory name
    cmd = [
        'cnmf', 'consensus',
        '--output-dir', str(dataset_dir.parent),
        '--name', dataset_dir.name,
        '--components', str(k_value),
        '--local-density-threshold', str(density_threshold),
        '--show-clustering'
    ]
    
    print(f"Command: {' '.join(cmd)}")
    
    try:
        # Run the consensus command
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print("STDOUT:", result.stdout)
        if result.stderr:
            print("STDERR:", result.stderr)
        
        # Format density threshold for filename (cNMF replaces '.' with '_')
        # Example: 0.01 -> 0_01, 0.1 -> 0_1
        dt_formatted = format_density_threshold(density_threshold)
        
        # The consensus file should be in the dataset directory
        # Format: {name}.consensus.k_{k}.dt_{dt_formatted}.consensus.txt
        consensus_file = dataset_dir / f"{dataset_dir.name}.consensus.k_{k_value}.dt_{dt_formatted}.consensus.txt"
        
        # The spectra file contains the GEP x gene matrix we need
        # cNMF can output with different naming conventions, try both:
        # 1. gene_spectra_score (Z-score normalized)
        # 2. gene_spectra_tpm (TPM normalized)
        spectra_file = dataset_dir / f"{dataset_dir.name}.gene_spectra_score.k_{k_value}.dt_{dt_formatted}.txt"
        
        if not spectra_file.exists():
            # Try TPM version
            spectra_file = dataset_dir / f"{dataset_dir.name}.gene_spectra_tpm.k_{k_value}.dt_{dt_formatted}.txt"
        
        if spectra_file.exists():
            print(f"✓ Consensus completed. Spectra file: {spectra_file.name}")
            return spectra_file
        else:
            # List available files to help debug
            print(f"\nERROR: Spectra file not found. Looking for:")
            print(f"  {dataset_dir.name}.gene_spectra_score.k_{k_value}.dt_{dt_formatted}.txt")
            print(f"  OR")
            print(f"  {dataset_dir.name}.gene_spectra_tpm.k_{k_value}.dt_{dt_formatted}.txt")
            print(f"\nAvailable files in {dataset_dir.name}:")
            for f in sorted(dataset_dir.glob(f"*k_{k_value}*.txt")):
                print(f"  {f.name}")
            raise FileNotFoundError(f"Spectra file not found in {dataset_dir}")
            
    except subprocess.CalledProcessError as e:
        print(f"ERROR running consensus: {e}")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        raise


def load_consensus_matrix(spectra_file, dataset_name):
    """
    Load the consensus spectra matrix (GEP x gene).
    
    Args:
        spectra_file (Path): Path to the spectra file
        dataset_name (str): Name of the dataset (for labeling)
        
    Returns:
        pd.DataFrame: Matrix with rows=GEPs, columns=genes
                     Index format: {dataset_name}_GEP_{i}
                     
    Purpose: Load the gene expression program matrix where each row represents
    a GEP (characterized by gene weights) and columns are genes. We rename
    the rows to include dataset name for tracking which patient each GEP comes from.
    """
    print(f"\nLoading consensus matrix from: {spectra_file.name}")
    
    # Load the spectra file - typically tab-separated with genes as columns
    df = pd.read_csv(spectra_file, sep='\t', index_col=0)
    
    # Rows should be GEPs (k programs), columns should be genes
    print(f"  Shape: {df.shape} (GEPs x genes)")
    print(f"  GEPs: {df.shape[0]}, Genes: {df.shape[1]}")
    
    # Rename index to include dataset name
    # Original index might be 0, 1, 2, ... or "GEP_0", "GEP_1", ...
    # We want: "patient1_GEP_0", "patient1_GEP_1", etc.
    df.index = [f"{dataset_name}_GEP_{i}" for i in range(len(df))]
    
    print(f"  Renamed GEPs: {df.index[0]} to {df.index[-1]}")
    
    return df


def stack_consensus_matrices(parent_dir, datasets, k_values, density_threshold):
    """
    Run consensus for all datasets and stack their GEP matrices vertically.
    
    Args:
        parent_dir (Path): Parent directory path
        datasets (list): List of dataset names
        k_values (list): List of k values (one per dataset)
        density_threshold (float): Density threshold for consensus
        
    Returns:
        pd.DataFrame: Stacked matrix where each row is a GEP from a dataset
        dict: Metadata about datasets (k values, number of GEPs)
        
    Purpose: Create a unified matrix containing all GEPs from all patients.
    Each row represents one GEP from one patient, enabling cross-patient
    comparison through clustering. Genes must be matched across datasets.
    """
    all_matrices = []
    metadata = {}
    
    for dataset, k_val in zip(datasets, k_values):
        dataset_dir = parent_dir / dataset
        
        # Run cNMF consensus for this dataset
        spectra_file = run_cnmf_consensus(dataset_dir, k_val, density_threshold)
        
        # Load the consensus matrix
        consensus_matrix = load_consensus_matrix(spectra_file, dataset)
        
        all_matrices.append(consensus_matrix)
        metadata[dataset] = {
            'k': k_val,
            'n_geps': consensus_matrix.shape[0],
            'n_genes': consensus_matrix.shape[1]
        }
    
    print(f"\n{'='*70}")
    print("Stacking consensus matrices...")
    print(f"{'='*70}")
    
    # Find common genes across all datasets
    # This is critical because different datasets may have different gene sets
    common_genes = set(all_matrices[0].columns)
    for matrix in all_matrices[1:]:
        common_genes = common_genes.intersection(set(matrix.columns))
    
    common_genes = sorted(list(common_genes))
    print(f"Common genes across all datasets: {len(common_genes)}")
    
    if len(common_genes) == 0:
        raise ValueError("No common genes found across datasets!")
    
    # Filter each matrix to common genes and stack vertically
    filtered_matrices = [matrix[common_genes] for matrix in all_matrices]
    stacked_matrix = pd.concat(filtered_matrices, axis=0)
    
    print(f"Stacked matrix shape: {stacked_matrix.shape}")
    print(f"  Total GEPs: {stacked_matrix.shape[0]}")
    print(f"  Common genes: {stacked_matrix.shape[1]}")
    print(f"\nGEP distribution by dataset:")
    for dataset, meta in metadata.items():
        print(f"  {dataset}: {meta['n_geps']} GEPs (k={meta['k']})")
    
    return stacked_matrix, metadata


def perform_hierarchical_clustering(stacked_matrix, linkage_method='average', 
                                    distance_metric='correlation'):
    """
    Perform agglomerative hierarchical clustering on GEPs.
    
    Args:
        stacked_matrix (pd.DataFrame): Matrix with rows=GEPs, columns=genes
        linkage_method (str): Linkage method ('average', 'complete', 'ward', etc.)
        distance_metric (str): Distance metric ('correlation', 'euclidean', etc.)
        
    Returns:
        np.ndarray: Linkage matrix from hierarchical clustering
        np.ndarray: Distance matrix between all GEPs
        
    Purpose: Cluster GEPs based on their gene expression patterns.
    - Distance metric: How to measure similarity between GEPs
      * 'correlation': 1 - Pearson correlation (common for expression data)
      * 'euclidean': Standard Euclidean distance
      * 'cosine': Cosine distance
    - Linkage method: How to merge clusters
      * 'average': UPGMA, good balance (recommended for expression data)
      * 'complete': Maximum distance (can create tight clusters)
      * 'ward': Minimizes within-cluster variance (good for clear groups)
    """
    print(f"\n{'='*70}")
    print("Performing Hierarchical Clustering")
    print(f"{'='*70}")
    print(f"Linkage method: {linkage_method}")
    print(f"Distance metric: {distance_metric}")
    
    # Calculate pairwise distances between GEPs
    # Each GEP is a row vector of gene weights
    print("\nCalculating pairwise distances...")
    distances = pdist(stacked_matrix.values, metric=distance_metric)
    
    # Convert to square distance matrix for visualization
    distance_matrix = squareform(distances)
    
    print(f"Distance matrix shape: {distance_matrix.shape}")
    print(f"Distance range: [{distance_matrix.min():.4f}, {distance_matrix.max():.4f}]")
    
    # Perform hierarchical clustering
    print(f"\nPerforming {linkage_method} linkage clustering...")
    linkage_matrix = linkage(distances, method=linkage_method)
    
    print(f"Linkage matrix shape: {linkage_matrix.shape}")
    print("✓ Clustering complete")
    
    return linkage_matrix, distance_matrix


def plot_dendrogram(linkage_matrix, labels, output_file, figsize=(20, 10)):
    """
    Create and save a dendrogram visualization.
    
    Args:
        linkage_matrix (np.ndarray): Linkage matrix from hierarchical clustering
        labels (list): GEP labels (dataset_GEP_i format)
        output_file (Path): Path to save the figure
        figsize (tuple): Figure size
        
    Purpose: Visualize the hierarchical relationships between GEPs.
    The dendrogram shows:
    - Height: Distance at which clusters merge
    - Leaves: Individual GEPs colored by dataset
    - Branches: Hierarchical grouping structure
    Helps identify which GEPs cluster together across patients.
    """
    print(f"\nCreating dendrogram...")
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Extract dataset names for coloring
    dataset_names = [label.split('_GEP_')[0] for label in labels]
    unique_datasets = sorted(set(dataset_names))
    
    # Create color map for datasets
    colors = sns.color_palette('husl', len(unique_datasets))
    dataset_colors = {dataset: colors[i] for i, dataset in enumerate(unique_datasets)}
    
    # Create dendrogram
    dendro = dendrogram(
        linkage_matrix,
        labels=labels,
        ax=ax,
        leaf_rotation=90,
        leaf_font_size=8,
        color_threshold=0
    )
    
    # Color the labels by dataset
    xlbls = ax.get_xmajorticklabels()
    for lbl in xlbls:
        dataset = lbl.get_text().split('_GEP_')[0]
        lbl.set_color(dataset_colors[dataset])
    
    ax.set_xlabel('Gene Expression Programs (GEPs)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Distance', fontsize=12, fontweight='bold')
    ax.set_title('Hierarchical Clustering of GEPs Across Datasets/Patients', 
                 fontsize=14, fontweight='bold', pad=20)
    
    # Add legend for datasets
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=dataset_colors[ds], label=ds) 
                      for ds in unique_datasets]
    ax.legend(handles=legend_elements, title='Dataset', 
             bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Dendrogram saved: {output_file}")
    plt.close()


def plot_heatmap(stacked_matrix, linkage_matrix, distance_matrix, output_file):
    """
    Create a clustered heatmap showing GEP similarities.
    
    Args:
        stacked_matrix (pd.DataFrame): Matrix with rows=GEPs, columns=genes
        linkage_matrix (np.ndarray): Linkage matrix for ordering
        distance_matrix (np.ndarray): Pairwise distance matrix
        output_file (Path): Path to save the figure
        
    Purpose: Visualize the distance/similarity matrix between all GEPs
    with hierarchical clustering-based ordering. This shows:
    - Blocks of similar GEPs (dark colors = low distance = high similarity)
    - Whether GEPs from the same patient cluster together (batch effects)
    - Cross-patient shared GEPs (off-diagonal blocks)
    """
    print(f"\nCreating clustered heatmap...")
    
    # Create figure with dendrogram
    fig = plt.figure(figsize=(16, 14))
    
    # Extract dataset names for row colors
    dataset_names = [label.split('_GEP_')[0] for label in stacked_matrix.index]
    unique_datasets = sorted(set(dataset_names))
    
    # Create color map for datasets
    colors = sns.color_palette('husl', len(unique_datasets))
    dataset_colors = {dataset: colors[i] for i, dataset in enumerate(unique_datasets)}
    row_colors = [dataset_colors[ds] for ds in dataset_names]
    
    # Create clustermap
    # This automatically reorders rows/columns based on hierarchical clustering
    g = sns.clustermap(
        distance_matrix,
        row_linkage=linkage_matrix,
        col_linkage=linkage_matrix,
        row_colors=row_colors,
        col_colors=row_colors,
        cmap='viridis',
        xticklabels=stacked_matrix.index,
        yticklabels=stacked_matrix.index,
        figsize=(16, 14),
        cbar_pos=(0.02, 0.8, 0.03, 0.15),
        dendrogram_ratio=0.15,
        cbar_kws={'label': 'Distance'}
    )
    
    # Adjust label size
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize=6, rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=6)
    
    # Add title
    g.fig.suptitle('GEP Distance Matrix with Hierarchical Clustering', 
                   fontsize=14, fontweight='bold', y=0.98)
    
    # Add dataset legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=dataset_colors[ds], label=ds) 
                      for ds in unique_datasets]
    g.ax_heatmap.legend(handles=legend_elements, title='Dataset',
                       bbox_to_anchor=(1.3, 1), loc='upper left', fontsize=9)
    
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Heatmap saved: {output_file}")
    plt.close()


def identify_gep_clusters(linkage_matrix, stacked_matrix, output_file, 
                         n_clusters=None, distance_threshold=None):
    """
    Cut the dendrogram to identify GEP clusters and save results.
    
    Args:
        linkage_matrix (np.ndarray): Linkage matrix
        stacked_matrix (pd.DataFrame): Original stacked matrix
        output_file (Path): Path to save cluster assignments
        n_clusters (int): Number of clusters to create (optional)
        distance_threshold (float): Distance threshold for cutting (optional)
        
    Returns:
        pd.DataFrame: Cluster assignments for each GEP
        
    Purpose: Group GEPs into discrete clusters to identify:
    - Shared programs: GEPs from different patients in the same cluster
    - Patient-specific programs: Clusters dominated by one patient
    - Cluster characteristics: Which genes define each cluster
    """
    print(f"\n{'='*70}")
    print("Identifying GEP Clusters")
    print(f"{'='*70}")
    
    # If neither n_clusters nor distance_threshold specified, use automatic cutting
    if n_clusters is None and distance_threshold is None:
        # Use the elbow method: cut at 70% of max distance
        max_distance = linkage_matrix[:, 2].max()
        distance_threshold = 0.7 * max_distance
        print(f"Using automatic distance threshold: {distance_threshold:.4f}")
    
    # Cut the dendrogram
    if distance_threshold is not None:
        clusters = fcluster(linkage_matrix, distance_threshold, criterion='distance')
    else:
        clusters = fcluster(linkage_matrix, n_clusters, criterion='maxclust')
    
    # Create results dataframe
    cluster_df = pd.DataFrame({
        'GEP': stacked_matrix.index,
        'Dataset': [label.split('_GEP_')[0] for label in stacked_matrix.index],
        'Cluster': clusters
    })
    
    print(f"Identified {cluster_df['Cluster'].nunique()} clusters")
    print("\nCluster summary:")
    print(cluster_df.groupby('Cluster')['Dataset'].value_counts().unstack(fill_value=0))
    
    # Save to file
    cluster_df.to_csv(output_file, index=False)
    print(f"\n✓ Cluster assignments saved: {output_file}")
    
    return cluster_df


def save_summary_statistics(stacked_matrix, linkage_matrix, metadata, 
                           cluster_df, output_file):
    """
    Save comprehensive summary statistics about the analysis.
    
    Args:
        stacked_matrix (pd.DataFrame): Stacked GEP matrix
        linkage_matrix (np.ndarray): Linkage matrix
        metadata (dict): Dataset metadata
        cluster_df (pd.DataFrame): Cluster assignments
        output_file (Path): Path to save summary
        
    Purpose: Create a human-readable summary of the analysis including:
    - Dataset information (samples, k values, GEP counts)
    - Clustering parameters and results
    - Cross-patient GEP sharing statistics
    - Cluster composition and characteristics
    """
    print(f"\nSaving summary statistics...")
    
    with open(output_file, 'w') as f:
        f.write("="*70 + "\n")
        f.write("cNMF CONSENSUS AND HIERARCHICAL CLUSTERING SUMMARY\n")
        f.write("="*70 + "\n\n")
        
        # Dataset information
        f.write("DATASET INFORMATION:\n")
        f.write("-" * 70 + "\n")
        for dataset, meta in metadata.items():
            f.write(f"{dataset}:\n")
            f.write(f"  K value: {meta['k']}\n")
            f.write(f"  Number of GEPs: {meta['n_geps']}\n")
            f.write(f"  Number of genes: {meta['n_genes']}\n\n")
        
        # Overall statistics
        f.write("\nOVERALL STATISTICS:\n")
        f.write("-" * 70 + "\n")
        f.write(f"Total number of GEPs: {stacked_matrix.shape[0]}\n")
        f.write(f"Common genes across datasets: {stacked_matrix.shape[1]}\n")
        f.write(f"Number of datasets: {len(metadata)}\n\n")
        
        # Clustering information
        f.write("\nCLUSTERING RESULTS:\n")
        f.write("-" * 70 + "\n")
        f.write(f"Number of clusters: {cluster_df['Cluster'].nunique()}\n")
        f.write(f"Clustering method: hierarchical agglomerative\n\n")
        
        f.write("Cluster composition (GEPs per dataset per cluster):\n")
        cluster_composition = cluster_df.groupby('Cluster')['Dataset'].value_counts().unstack(fill_value=0)
        f.write(cluster_composition.to_string())
        f.write("\n\n")
        
        # Shared vs patient-specific clusters
        f.write("\nSHARED vs PATIENT-SPECIFIC ANALYSIS:\n")
        f.write("-" * 70 + "\n")
        for cluster_id in sorted(cluster_df['Cluster'].unique()):
            cluster_data = cluster_df[cluster_df['Cluster'] == cluster_id]
            n_datasets = cluster_data['Dataset'].nunique()
            total_geps = len(cluster_data)
            
            if n_datasets > 1:
                cluster_type = "SHARED across patients"
            else:
                cluster_type = f"PATIENT-SPECIFIC ({cluster_data['Dataset'].iloc[0]})"
            
            f.write(f"\nCluster {cluster_id} ({cluster_type}):\n")
            f.write(f"  Total GEPs: {total_geps}\n")
            f.write(f"  Datasets involved: {n_datasets}\n")
            f.write(f"  Dataset distribution: {dict(cluster_data['Dataset'].value_counts())}\n")
    
    print(f"✓ Summary saved: {output_file}")


def main():
    """
    Main execution function orchestrating the entire analysis workflow.
    
    Workflow:
    1. Parse command line arguments
    2. Identify dataset directories
    3. Run cNMF consensus for each dataset
    4. Load and stack consensus matrices
    5. Perform hierarchical clustering
    6. Create visualizations
    7. Identify and characterize clusters
    8. Save results and summary
    """
    # Parse arguments
    args = parse_arguments()
    
    print("\n" + "="*70)
    print("cNMF CONSENSUS AND HIERARCHICAL CLUSTERING ANALYSIS")
    print("="*70)
    print(f"\nConfiguration:")
    print(f"  Parent directory: {args.parent_dir}")
    print(f"  Selected k values: {args.selected_ks}")
    print(f"  Density threshold: {args.density_threshold}")
    print(f"  Output prefix: {args.output_prefix}")
    print(f"  Linkage method: {args.linkage_method}")
    print(f"  Distance metric: {args.distance_metric}")
    
    # Set up paths
    parent_dir = Path(args.parent_dir)
    
    # Parse k values
    k_values = [int(k.strip()) for k in args.selected_ks.split(',')]
    
    # Get dataset directories
    datasets = get_dataset_directories(parent_dir)
    
    # Validate k values match datasets
    if len(k_values) != len(datasets):
        raise ValueError(f"Number of k values ({len(k_values)}) must match "
                        f"number of datasets ({len(datasets)})")
    
    # Stack consensus matrices from all datasets
    stacked_matrix, metadata = stack_consensus_matrices(
        parent_dir, datasets, k_values, args.density_threshold
    )
    
    # Save the stacked matrix
    stacked_matrix_file = f"{args.output_prefix}_stacked_gep_matrix.csv"
    stacked_matrix.to_csv(stacked_matrix_file)
    print(f"\n✓ Stacked GEP matrix saved: {stacked_matrix_file}")
    
    # Perform hierarchical clustering
    linkage_matrix, distance_matrix = perform_hierarchical_clustering(
        stacked_matrix, 
        linkage_method=args.linkage_method,
        distance_metric=args.distance_metric
    )
    
    # Create visualizations
    dendrogram_file = f"{args.output_prefix}_dendrogram.png"
    plot_dendrogram(linkage_matrix, stacked_matrix.index.tolist(), 
                   dendrogram_file, figsize=(24, 12))
    
    heatmap_file = f"{args.output_prefix}_distance_heatmap.png"
    plot_heatmap(stacked_matrix, linkage_matrix, distance_matrix, heatmap_file)
    
    # Identify clusters
    cluster_file = f"{args.output_prefix}_cluster_assignments.csv"
    cluster_df = identify_gep_clusters(linkage_matrix, stacked_matrix, cluster_file)
    
    # Save summary statistics
    summary_file = f"{args.output_prefix}_analysis_summary.txt"
    save_summary_statistics(stacked_matrix, linkage_matrix, metadata, 
                           cluster_df, summary_file)
    
    # Final summary
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE!")
    print("="*70)
    print("\nOutput files generated:")
    print(f"  1. {stacked_matrix_file} - Stacked GEP x gene matrix")
    print(f"  2. {dendrogram_file} - Hierarchical clustering dendrogram")
    print(f"  3. {heatmap_file} - Clustered distance heatmap")
    print(f"  4. {cluster_file} - Cluster assignments for each GEP")
    print(f"  5. {summary_file} - Comprehensive analysis summary")
    print("\nInterpretation guide:")
    print("  • Dendrogram: Shows hierarchical relationships between GEPs")
    print("  • Heatmap: Dark regions = similar GEPs, off-diagonal blocks = cross-patient similarity")
    print("  • Clusters: Groups of similar GEPs that may represent shared programs")
    print("  • Summary: Statistics about patient-specific vs shared GEPs")
    print("="*70 + "\n")


if __name__ == '__main__':
    main()