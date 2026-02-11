#!/usr/bin/env python3
"""
cNMF K-Value Comparison: Hierarchical Clustering Analysis
==========================================================

This script analyzes gene expression programs (GEPs) identified by cNMF across
different k values for a SINGLE dataset to identify:
- Which GEPs are consistently found across different k values
- How GEPs cluster together across k values
- Correlation structure of GEPs across k values

Author: Analysis pipeline for cNMF consensus comparison
"""

import argparse
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from scipy.stats import zscore
import warnings
import cnmf
warnings.filterwarnings('ignore')


def run_consensus_for_k(cnmf_dir, k_value, density_threshold=0.01):
    """
    Run cNMF consensus for a specific k value.
    
    Parameters:
    -----------
    cnmf_dir : str or Path
        Path to the cNMF results directory
    k_value : int
        Number of components (k)
    density_threshold : float
        Density threshold for consensus (default: 0.01)
    
    Returns:
    --------
    str : Path to the gene spectra score file
    """
    print(f"  Running consensus for k={k_value}...")
    
    # Initialize cNMF object
    cnmf_obj = cnmf.cNMF(output_dir=str(cnmf_dir.parent), name=Path(cnmf_dir).name)
    
    # Run consensus
    try:
        cnmf_obj.consensus(k=k_value, 
                          density_threshold=density_threshold,
                          show_clustering=False, refit_usage=False)
        
        # Construct expected output file path
        # NOTE: cNMF replaces decimal point with underscore in filenames
        # e.g., dt_0.01 becomes dt_0_01
        dt_str = f"{density_threshold:.2f}".replace('.', '_')
        gene_spectra_file = os.path.join(
            cnmf_dir,
            f"{Path(cnmf_dir).name}.gene_spectra_score.k_{k_value}.dt_{dt_str}.txt"
        )
        
        if os.path.exists(gene_spectra_file):
            print(f"    ✓ Consensus completed. Output: {os.path.basename(gene_spectra_file)}")
            return gene_spectra_file
        else:
            print(f"    ✗ Expected output file not found: {gene_spectra_file}")
            return None
            
    except Exception as e:
        print(f"    ✗ Error running consensus for k={k_value}: {str(e)}")
        return None


def load_gene_spectra(file_path, k_value, dataset_name):
    """
    Load gene spectra score matrix and add metadata.
    
    Parameters:
    -----------
    file_path : str
        Path to gene spectra score file
    k_value : int
        K value for this matrix
    dataset_name : str
        Name of the dataset
    
    Returns:
    --------
    pd.DataFrame : Gene spectra matrix with renamed indices
    """
    print(f"  Loading gene spectra for k={k_value}...")
    
    try:
        # Load the matrix (rows = GEPs, columns = genes)
        df = pd.read_csv(file_path, sep='\t', index_col=0)
        
        # Rename indices to include k value: k10_GEP1, k10_GEP2, etc.
        df.index = [f"k{k_value}_GEP{i+1}" for i in range(len(df))]
        
        print(f"    ✓ Loaded {len(df)} GEPs × {len(df.columns)} genes")
        return df
        
    except Exception as e:
        print(f"    ✗ Error loading file: {str(e)}")
        return None


def stack_matrices(gene_spectra_dict):
    """
    Stack all gene spectra matrices vertically.
    
    Parameters:
    -----------
    gene_spectra_dict : dict
        Dictionary mapping k values to gene spectra DataFrames
    
    Returns:
    --------
    pd.DataFrame : Stacked matrix with all GEPs from all k values
    """
    print("\nStacking gene spectra matrices...")
    
    # Get all k values in sorted order
    k_values = sorted(gene_spectra_dict.keys())
    
    # Stack matrices vertically
    stacked = pd.concat([gene_spectra_dict[k] for k in k_values], axis=0)
    
    print(f"  ✓ Stacked matrix: {len(stacked)} total GEPs × {len(stacked.columns)} genes")
    print(f"  ✓ K values included: {', '.join(map(str, k_values))}")
    
    # Print breakdown by k
    for k in k_values:
        n_geps = len(gene_spectra_dict[k])
        print(f"    - k={k}: {n_geps} GEPs")
    
    return stacked


def calculate_correlation_matrix(stacked_matrix):
    """
    Calculate Pearson correlation matrix between all GEPs.
    
    Parameters:
    -----------
    stacked_matrix : pd.DataFrame
        Stacked gene spectra matrix
    
    Returns:
    --------
    pd.DataFrame : Pearson correlation matrix (GEP × GEP)
    """
    print("\nCalculating Pearson correlation matrix...")
    
    # Calculate correlation between rows (GEPs)
    # Each row is a GEP, each column is a gene
    corr_matrix = stacked_matrix.T.corr(method='pearson')
    
    print(f"  ✓ Correlation matrix: {corr_matrix.shape[0]} × {corr_matrix.shape[1]}")
    
    # Print some statistics
    # Exclude diagonal (self-correlation = 1.0)
    off_diagonal = corr_matrix.values[~np.eye(corr_matrix.shape[0], dtype=bool)]
    print(f"  ✓ Correlation statistics (excluding diagonal):")
    print(f"    - Mean: {np.mean(off_diagonal):.3f}")
    print(f"    - Median: {np.median(off_diagonal):.3f}")
    print(f"    - Min: {np.min(off_diagonal):.3f}")
    print(f"    - Max: {np.max(off_diagonal):.3f}")
    
    return corr_matrix


def perform_hierarchical_clustering(correlation_matrix, linkage_method='average'):
    """
    Perform hierarchical clustering on correlation matrix.
    
    Parameters:
    -----------
    correlation_matrix : pd.DataFrame
        Pearson correlation matrix
    linkage_method : str
        Linkage method for clustering (default: 'average')
    
    Returns:
    --------
    np.ndarray : Linkage matrix from hierarchical clustering
    """
    print(f"\nPerforming hierarchical clustering (method: {linkage_method})...")
    
    # Convert correlation to distance: distance = 1 - correlation
    # Higher correlation = lower distance
    distance_matrix = 1 - correlation_matrix.values
    
    # Ensure distance matrix is valid (symmetric, non-negative, zero diagonal)
    distance_matrix = np.clip(distance_matrix, 0, 2)  # Correlation ranges from -1 to 1
    np.fill_diagonal(distance_matrix, 0)
    
    # Convert to condensed distance matrix for scipy
    condensed_dist = squareform(distance_matrix, checks=False)
    
    # Perform hierarchical clustering
    linkage_matrix = linkage(condensed_dist, method=linkage_method)
    
    print(f"  ✓ Hierarchical clustering completed")
    
    return linkage_matrix


def plot_dendrogram(linkage_matrix, labels, output_prefix, figsize=(14, 8)):
    """
    Create and save dendrogram visualization.
    
    Parameters:
    -----------
    linkage_matrix : np.ndarray
        Linkage matrix from hierarchical clustering
    labels : list
        Labels for each GEP (e.g., ['k10_GEP1', 'k10_GEP2', ...])
    output_prefix : str
        Prefix for output file
    figsize : tuple
        Figure size (width, height)
    """
    print("\nCreating dendrogram...")
    
    # Extract k values from labels for coloring
    k_values = [int(label.split('_')[0][1:]) for label in labels]
    unique_k = sorted(set(k_values))
    
    # Create color map
    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_k)))
    k_to_color = {k: colors[i] for i, k in enumerate(unique_k)}
    label_colors = [k_to_color[k] for k in k_values]
    
    # Create figure
    plt.figure(figsize=figsize)
    
    # Plot dendrogram
    dendro = dendrogram(
        linkage_matrix,
        labels=labels,
        leaf_rotation=90,
        leaf_font_size=8,
        color_threshold=0
    )
    
    # Color the labels by k value
    ax = plt.gca()
    xlbls = ax.get_xmajorticklabels()
    for lbl in xlbls:
        label_text = lbl.get_text()
        if label_text in labels:
            idx = labels.index(label_text)
            lbl.set_color(label_colors[idx])
    
    plt.title('Hierarchical Clustering of GEPs Across K Values', fontsize=14, fontweight='bold')
    plt.xlabel('GEP ID', fontsize=12)
    plt.ylabel('Distance (1 - Pearson Correlation)', fontsize=12)
    
    # Add legend
    legend_elements = [plt.Line2D([0], [0], color=k_to_color[k], lw=4, label=f'k={k}')
                      for k in unique_k]
    plt.legend(handles=legend_elements, loc='upper right', title='K Value')
    
    plt.tight_layout()
    
    output_file = f"{output_prefix}_dendrogram.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Dendrogram saved: {output_file}")
    plt.close()


def plot_correlation_heatmap(correlation_matrix, output_prefix, figsize=(12, 10)):
    """
    Create and save correlation heatmap with hierarchical clustering.
    Now includes color-coded axes showing k-value for each GEP.
    
    Parameters:
    -----------
    correlation_matrix : pd.DataFrame
        Pearson correlation matrix
    output_prefix : str
        Prefix for output file
    figsize : tuple
        Figure size (width, height)
    """
    print("\nCreating correlation heatmap...")
    
    # Extract k values from row/column labels
    labels = correlation_matrix.index.tolist()
    k_values = [int(label.split('_')[0][1:]) for label in labels]
    unique_k = sorted(set(k_values))
    
    # Create color map for k values
    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_k)))
    k_to_color = {k: colors[i] for i, k in enumerate(unique_k)}
    
    # Create row and column colors for the heatmap
    row_colors = [k_to_color[k] for k in k_values]
    col_colors = [k_to_color[k] for k in k_values]
    
    # Use seaborn clustermap for automatic clustering with color bars
    g = sns.clustermap(
        correlation_matrix,
        cmap='RdBu_r',
        center=0,
        vmin=-1, vmax=1,
        figsize=figsize,
        dendrogram_ratio=0.15,
        cbar_pos=(0.02, 0.83, 0.03, 0.15),
        linewidths=0,
        xticklabels=True,
        yticklabels=True,
        cbar_kws={'label': 'Pearson Correlation'},
        row_colors=row_colors,
        col_colors=col_colors
    )
    
    # Rotate labels
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), rotation=90, fontsize=6)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), rotation=0, fontsize=6)
    
    # Add title
    plt.suptitle('Pearson Correlation Between GEPs Across K Values', 
                 fontsize=14, fontweight='bold', y=0.98)
    
    # Add legend for k-value colors
    # Create legend patches
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=k_to_color[k], label=f'k={k}')
                      for k in unique_k]
    
    # Add legend to the right side of the plot
    g.ax_heatmap.legend(handles=legend_elements, 
                       title='K Value', 
                       bbox_to_anchor=(1.35, 1.0),
                       loc='upper left',
                       frameon=True)
    
    output_file = f"{output_prefix}_correlation_heatmap.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Correlation heatmap saved: {output_file}")
    plt.close()


def plot_pairwise_k_comparison(correlation_matrix, k1, k2, output_prefix):
    """
    Create pairwise comparison heatmap between GEPs from two different k values.
    
    Parameters:
    -----------
    correlation_matrix : pd.DataFrame
        Full Pearson correlation matrix
    k1 : int
        First k value
    k2 : int
        Second k value
    output_prefix : str
        Prefix for output file
    """
    print(f"  Creating pairwise comparison heatmap: k={k1} vs k={k2}...")
    
    # Get GEPs for each k value
    k1_geps = [gep for gep in correlation_matrix.index if gep.startswith(f'k{k1}_')]
    k2_geps = [gep for gep in correlation_matrix.index if gep.startswith(f'k{k2}_')]
    
    if not k1_geps or not k2_geps:
        print(f"    ✗ No GEPs found for k={k1} or k={k2}")
        return
    
    # Extract the submatrix for this pair
    pairwise_corr = correlation_matrix.loc[k1_geps, k2_geps]
    
    # Calculate figure size based on number of GEPs
    n_k1 = len(k1_geps)
    n_k2 = len(k2_geps)
    figsize = (max(8, n_k2 * 0.5), max(6, n_k1 * 0.5))
    
    # Create heatmap
    plt.figure(figsize=figsize)
    
    # Plot heatmap with better visibility
    sns.heatmap(
        pairwise_corr,
        cmap='RdBu_r',
        center=0,
        vmin=-1, vmax=1,
        annot=True if (n_k1 <= 12 and n_k2 <= 12) else False,  # Only annotate if small
        fmt='.2f',
        cbar_kws={'label': 'Pearson Correlation'},
        linewidths=0.5,
        linecolor='gray',
        square=False
    )
    
    plt.title(f'Pairwise GEP Correlation: k={k1} vs k={k2}', 
             fontsize=14, fontweight='bold', pad=20)
    plt.xlabel(f'GEPs from k={k2}', fontsize=12, fontweight='bold')
    plt.ylabel(f'GEPs from k={k1}', fontsize=12, fontweight='bold')
    plt.tight_layout()
    
    output_file = f"{output_prefix}_pairwise_k{k1}_vs_k{k2}.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"    ✓ Pairwise heatmap saved: {output_file}")
    plt.close()


def create_all_pairwise_comparisons(correlation_matrix, k_values, output_prefix):
    """
    Create pairwise comparison heatmaps for all combinations of k values.
    Only called when there are 4 or fewer k values.
    
    Parameters:
    -----------
    correlation_matrix : pd.DataFrame
        Full Pearson correlation matrix
    k_values : list
        List of k values to compare
    output_prefix : str
        Prefix for output files
    """
    print("\nCreating pairwise k-value comparison heatmaps...")
    print(f"  Number of k values: {len(k_values)}")
    
    k_sorted = sorted(k_values)
    n_comparisons = 0
    
    # Generate all pairwise combinations
    for i in range(len(k_sorted)):
        for j in range(i + 1, len(k_sorted)):
            k1 = k_sorted[i]
            k2 = k_sorted[j]
            plot_pairwise_k_comparison(correlation_matrix, k1, k2, output_prefix)
            n_comparisons += 1
    
    print(f"\n  ✓ Created {n_comparisons} pairwise comparison heatmaps")
    print(f"  ✓ Files: {output_prefix}_pairwise_k*_vs_k*.png")



def identify_similar_geps(correlation_matrix, threshold=0.8, output_prefix=None):
    """
    Identify highly correlated GEPs across different k values.
    
    Parameters:
    -----------
    correlation_matrix : pd.DataFrame
        Pearson correlation matrix
    threshold : float
        Correlation threshold for identifying similar GEPs (default: 0.8)
    output_prefix : str, optional
        If provided, save results to file
    
    Returns:
    --------
    pd.DataFrame : Table of highly correlated GEP pairs
    """
    print(f"\nIdentifying similar GEPs (correlation ≥ {threshold})...")
    
    similar_pairs = []
    
    # Get upper triangle indices (avoid duplicates and self-correlations)
    rows, cols = np.triu_indices_from(correlation_matrix, k=1)
    
    for i, j in zip(rows, cols):
        gep1 = correlation_matrix.index[i]
        gep2 = correlation_matrix.index[j]
        corr = correlation_matrix.iloc[i, j]
        
        if corr >= threshold:
            # Extract k values
            k1 = int(gep1.split('_')[0][1:])
            k2 = int(gep2.split('_')[0][1:])
            
            similar_pairs.append({
                'GEP_1': gep1,
                'GEP_2': gep2,
                'K_1': k1,
                'K_2': k2,
                'Correlation': corr,
                'Cross_K': k1 != k2
            })
    
    if similar_pairs:
        df = pd.DataFrame(similar_pairs)
        df = df.sort_values('Correlation', ascending=False)
        
        print(f"  ✓ Found {len(df)} GEP pairs with correlation ≥ {threshold}")
        
        # Count cross-k vs same-k pairs
        cross_k = df['Cross_K'].sum()
        same_k = len(df) - cross_k
        print(f"    - Cross-k pairs (GEPs from different k): {cross_k}")
        print(f"    - Same-k pairs (GEPs from same k): {same_k}")
        
        if output_prefix:
            output_file = f"{output_prefix}_similar_geps_threshold{threshold:.2f}.csv"
            df.to_csv(output_file, index=False)
            print(f"  ✓ Similar GEPs table saved: {output_file}")
        
        return df
    else:
        print(f"  ✗ No GEP pairs found with correlation ≥ {threshold}")
        return pd.DataFrame()


def create_summary_report(stacked_matrix, correlation_matrix, gene_spectra_dict, 
                         linkage_matrix, output_prefix, args):
    """
    Create a comprehensive text summary of the analysis.
    
    Parameters:
    -----------
    stacked_matrix : pd.DataFrame
        Stacked gene spectra matrix
    correlation_matrix : pd.DataFrame
        Pearson correlation matrix
    gene_spectra_dict : dict
        Dictionary of gene spectra matrices by k value
    linkage_matrix : np.ndarray
        Hierarchical clustering linkage matrix
    output_prefix : str
        Prefix for output file
    args : argparse.Namespace
        Command line arguments
    """
    print("\nGenerating summary report...")
    
    output_file = f"{output_prefix}_analysis_summary.txt"
    
    with open(output_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("cNMF K-VALUE COMPARISON: HIERARCHICAL CLUSTERING ANALYSIS\n")
        f.write("="*80 + "\n\n")
        
        # Analysis parameters
        f.write("ANALYSIS PARAMETERS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Dataset directory: {args.cnmf_dir}\n")
        f.write(f"Dataset name: {Path(args.cnmf_dir).name}\n")
        f.write(f"K values analyzed: {', '.join(map(str, sorted(gene_spectra_dict.keys())))}\n")
        f.write(f"Density threshold: {args.density_threshold}\n")
        f.write(f"Linkage method: {args.linkage_method}\n")
        f.write(f"Number of genes: {stacked_matrix.shape[1]}\n\n")
        
        # GEP counts
        f.write("GEP COUNTS BY K VALUE\n")
        f.write("-" * 80 + "\n")
        total_geps = 0
        for k in sorted(gene_spectra_dict.keys()):
            n_geps = len(gene_spectra_dict[k])
            total_geps += n_geps
            f.write(f"  k={k:3d}: {n_geps:3d} GEPs\n")
        f.write(f"  Total: {total_geps:3d} GEPs\n\n")
        
        # Correlation statistics
        f.write("CORRELATION STATISTICS\n")
        f.write("-" * 80 + "\n")
        off_diag = correlation_matrix.values[~np.eye(correlation_matrix.shape[0], dtype=bool)]
        f.write(f"  Mean correlation: {np.mean(off_diag):.4f}\n")
        f.write(f"  Median correlation: {np.median(off_diag):.4f}\n")
        f.write(f"  Std deviation: {np.std(off_diag):.4f}\n")
        f.write(f"  Min correlation: {np.min(off_diag):.4f}\n")
        f.write(f"  Max correlation: {np.max(off_diag):.4f}\n\n")
        
        # Correlation distribution
        f.write("CORRELATION DISTRIBUTION\n")
        f.write("-" * 80 + "\n")
        bins = [(-1.0, -0.5), (-0.5, 0.0), (0.0, 0.3), (0.3, 0.5), 
                (0.5, 0.7), (0.7, 0.9), (0.9, 1.0)]
        for low, high in bins:
            count = np.sum((off_diag >= low) & (off_diag < high))
            pct = 100 * count / len(off_diag)
            f.write(f"  [{low:5.2f}, {high:5.2f}): {count:6d} pairs ({pct:5.2f}%)\n")
        f.write("\n")
        
        # High correlation pairs analysis
        f.write("HIGHLY CORRELATED GEP PAIRS (correlation ≥ 0.8)\n")
        f.write("-" * 80 + "\n")
        
        similar_df = identify_similar_geps(correlation_matrix, threshold=0.8, 
                                          output_prefix=None)
        
        if len(similar_df) > 0:
            cross_k = similar_df['Cross_K'].sum()
            same_k = len(similar_df) - cross_k
            
            f.write(f"  Total highly correlated pairs: {len(similar_df)}\n")
            f.write(f"  Cross-k pairs (different k values): {cross_k}\n")
            f.write(f"  Same-k pairs (same k value): {same_k}\n\n")
            
            if cross_k > 0:
                f.write("  Cross-k highly correlated pairs (top 20):\n")
                cross_k_pairs = similar_df[similar_df['Cross_K']].head(20)
                for idx, row in cross_k_pairs.iterrows():
                    f.write(f"    {row['GEP_1']:15s} ↔ {row['GEP_2']:15s}  "
                           f"(r = {row['Correlation']:.4f})\n")
        else:
            f.write("  No highly correlated pairs found (threshold = 0.8)\n")
        f.write("\n")
        
        # Interpretation guide
        f.write("="*80 + "\n")
        f.write("INTERPRETATION GUIDE\n")
        f.write("="*80 + "\n\n")
        
        f.write("WHAT DO THESE RESULTS MEAN?\n")
        f.write("-" * 80 + "\n")
        f.write("This analysis helps you understand:\n\n")
        
        f.write("1. CONSISTENCY ACROSS K VALUES:\n")
        f.write("   • High cross-k correlations indicate GEPs that are consistently\n")
        f.write("     identified regardless of the k parameter choice\n")
        f.write("   • These are likely robust, biologically meaningful programs\n\n")
        
        f.write("2. K VALUE STABILITY:\n")
        f.write("   • If GEPs from adjacent k values (e.g., k=10 and k=11) cluster\n")
        f.write("     together, it suggests those k values capture similar biology\n")
        f.write("   • Large jumps in GEP composition between k values may indicate\n")
        f.write("     a transition in model complexity\n\n")
        
        f.write("3. OPTIMAL K SELECTION:\n")
        f.write("   • Look for k values where:\n")
        f.write("     - GEPs are distinct from each other (low intra-k correlation)\n")
        f.write("     - But robust programs are preserved across k (high cross-k correlation)\n")
        f.write("   • The dendrogram can help visualize when new GEPs emerge vs when\n")
        f.write("     existing GEPs split into sub-programs\n\n")
        
        f.write("HOW TO USE THE OUTPUT FILES:\n")
        f.write("-" * 80 + "\n")
        f.write("• Dendrogram: Shows hierarchical relationships between ALL GEPs\n")
        f.write("  - GEPs are colored by k value\n")
        f.write("  - Mixed-color clusters = programs stable across k values\n")
        f.write("  - Single-color clusters = k-specific programs\n\n")
        
        f.write("• Correlation Heatmap: Shows pairwise correlations\n")
        f.write("  - Dark red blocks = highly correlated GEPs\n")
        f.write("  - Off-diagonal blocks = cross-k correlations\n")
        f.write("  - Use to identify which specific GEPs are similar\n\n")
        
        f.write("• Similar GEPs Table: Lists all high-correlation pairs\n")
        f.write("  - Focus on cross-k pairs to find stable programs\n")
        f.write("  - These GEPs likely represent the same biological process\n\n")
        
        f.write("• Stacked Matrix: Full gene weights for all GEPs\n")
        f.write("  - Use for downstream analysis (e.g., gene set enrichment)\n")
        f.write("  - Compare gene loadings for similar GEPs across k values\n\n")
        
        f.write("="*80 + "\n")
        f.write("NEXT STEPS\n")
        f.write("="*80 + "\n")
        f.write("1. Examine the dendrogram and heatmap to identify stable GEPs\n")
        f.write("2. For highly correlated cross-k GEPs, compare their top genes\n")
        f.write("3. Consider the k value that:\n")
        f.write("   a) Captures the key stable programs (seen across multiple k)\n")
        f.write("   b) Provides useful resolution of biological processes\n")
        f.write("   c) Doesn't over-split programs into redundant components\n")
        f.write("4. Use biological knowledge and gene set enrichment to validate\n")
        f.write("   that the selected k captures meaningful programs\n")
        f.write("="*80 + "\n")
    
    print(f"  ✓ Summary report saved: {output_file}")


def main():
    """Main analysis pipeline."""
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Hierarchical clustering analysis of cNMF GEPs across k values',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python cnmf_k_comparison_clustering.py \\
      --cnmf_dir /path/to/dataset \\
      --k_values 8 9 10 11 12 \\
      --density_threshold 0.01 \\
      --output_prefix dataset_k_comparison

This will:
  1. Run cNMF consensus for k=8,9,10,11,12
  2. Load gene spectra score matrices for each k
  3. Stack all matrices vertically
  4. Calculate Pearson correlation between all GEPs
  5. Perform hierarchical clustering
  6. Generate visualizations and summary report
        """
    )
    
    parser.add_argument('--cnmf_dir', required=True,
                       help='Path to cNMF results directory for the dataset')
    parser.add_argument('--k_values', required=True, nargs='+', type=int,
                       help='K values to analyze (space-separated)')
    parser.add_argument('--density_threshold', type=float, default=0.01,
                       help='Density threshold for cNMF consensus (default: 0.01)')
    parser.add_argument('--linkage_method', default='average',
                       choices=['average', 'complete', 'single', 'ward', 'weighted'],
                       help='Linkage method for hierarchical clustering (default: average)')
    parser.add_argument('--output_prefix', default='k_comparison',
                       help='Prefix for output files (default: k_comparison)')
    parser.add_argument('--skip_consensus', action='store_true',
                       help='Skip running consensus (assumes already completed)')
    
    args = parser.parse_args()
    
    # Validate inputs
    cnmf_dir = Path(args.cnmf_dir)
    if not cnmf_dir.exists():
        print(f"ERROR: Directory not found: {cnmf_dir}")
        sys.exit(1)
    
    if len(args.k_values) == 0:
        print("ERROR: Must provide at least one k value")
        sys.exit(1)
    
    # Print header
    print("="*80)
    print("cNMF K-VALUE COMPARISON: HIERARCHICAL CLUSTERING ANALYSIS")
    print("="*80)
    print(f"\nDataset: {cnmf_dir.name}")
    print(f"K values: {', '.join(map(str, sorted(args.k_values)))}")
    print(f"Density threshold: {args.density_threshold}")
    print(f"Linkage method: {args.linkage_method}")
    print("\n" + "="*80)
    
    # Step 1: Run consensus for each k value (unless skipped)
    gene_spectra_files = {}
    
    if not args.skip_consensus:
        print("\nSTEP 1: Running cNMF consensus for each k value")
        print("-"*80)
        
        for k in sorted(args.k_values):
            gene_spectra_file = run_consensus_for_k(
                cnmf_dir, k, args.density_threshold
            )
            if gene_spectra_file:
                gene_spectra_files[k] = gene_spectra_file
            else:
                print(f"  WARNING: Consensus failed for k={k}, skipping...")
    else:
        print("\nSTEP 1: Skipping consensus (--skip_consensus flag set)")
        print("-"*80)
        print("Looking for existing gene spectra files...")
        
        dataset_name = cnmf_dir.name
        # NOTE: cNMF replaces decimal point with underscore in filenames
        dt_str = f"{args.density_threshold:.2f}".replace('.', '_')
        
        for k in sorted(args.k_values):
            expected_file = cnmf_dir / f"{dataset_name}.gene_spectra_score.k_{k}.dt_{dt_str}.txt"
            if expected_file.exists():
                gene_spectra_files[k] = str(expected_file)
                print(f"  ✓ Found k={k}: {expected_file.name}")
            else:
                print(f"  ✗ Not found k={k}: {expected_file.name}")
    
    if not gene_spectra_files:
        print("\nERROR: No gene spectra files found. Cannot proceed.")
        sys.exit(1)
    
    # Step 2: Load gene spectra matrices
    print("\n" + "="*80)
    print("STEP 2: Loading gene spectra matrices")
    print("-"*80)
    
    gene_spectra_dict = {}
    dataset_name = cnmf_dir.name
    
    for k, file_path in sorted(gene_spectra_files.items()):
        df = load_gene_spectra(file_path, k, dataset_name)
        if df is not None:
            gene_spectra_dict[k] = df
    
    if not gene_spectra_dict:
        print("\nERROR: Failed to load any gene spectra matrices. Cannot proceed.")
        sys.exit(1)
    
    # Step 3: Stack matrices
    print("\n" + "="*80)
    print("STEP 3: Stacking gene spectra matrices")
    print("-"*80)
    
    stacked_matrix = stack_matrices(gene_spectra_dict)
    
    # Save stacked matrix
    output_file = f"{args.output_prefix}_stacked_gep_matrix.csv"
    stacked_matrix.to_csv(output_file)
    print(f"\n  ✓ Stacked matrix saved: {output_file}")
    
    # Step 4: Calculate correlation matrix
    print("\n" + "="*80)
    print("STEP 4: Calculating Pearson correlation matrix")
    print("-"*80)
    
    correlation_matrix = calculate_correlation_matrix(stacked_matrix)
    
    # Save correlation matrix
    output_file = f"{args.output_prefix}_correlation_matrix.csv"
    correlation_matrix.to_csv(output_file)
    print(f"\n  ✓ Correlation matrix saved: {output_file}")
    
    # Step 5: Hierarchical clustering
    print("\n" + "="*80)
    print("STEP 5: Performing hierarchical clustering")
    print("-"*80)
    
    linkage_matrix = perform_hierarchical_clustering(
        correlation_matrix, 
        args.linkage_method
    )
    
    # Step 6: Create visualizations
    print("\n" + "="*80)
    print("STEP 6: Creating visualizations")
    print("-"*80)
    
    plot_dendrogram(
        linkage_matrix,
        list(stacked_matrix.index),
        args.output_prefix
    )
    
    plot_correlation_heatmap(
        correlation_matrix,
        args.output_prefix
    )
    
    # Create pairwise comparison heatmaps if 4 or fewer k values
    if len(gene_spectra_dict) <= 4:
        print("\n" + "-"*80)
        print("Creating pairwise k-value comparison heatmaps...")
        print(f"(Generated because number of k values [{len(gene_spectra_dict)}] ≤ 4)")
        print("-"*80)
        create_all_pairwise_comparisons(
            correlation_matrix,
            list(gene_spectra_dict.keys()),
            args.output_prefix
        )
    else:
        print("\n" + "-"*80)
        print(f"Skipping pairwise comparison heatmaps (number of k values [{len(gene_spectra_dict)}] > 4)")
        print("These are only generated when comparing 4 or fewer k values")
        print("-"*80)
    
    # Step 7: Identify similar GEPs
    print("\n" + "="*80)
    print("STEP 7: Identifying similar GEPs")
    print("-"*80)
    
    for threshold in [0.9, 0.8, 0.7]:
        identify_similar_geps(
            correlation_matrix,
            threshold=threshold,
            output_prefix=args.output_prefix
        )
    
    # Step 8: Create summary report
    print("\n" + "="*80)
    print("STEP 8: Creating summary report")
    print("-"*80)
    
    create_summary_report(
        stacked_matrix,
        correlation_matrix,
        gene_spectra_dict,
        linkage_matrix,
        args.output_prefix,
        args
    )
    
    # Final summary
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80)
    print("\nOutput files generated:")
    print(f"  1. {args.output_prefix}_stacked_gep_matrix.csv")
    print(f"  2. {args.output_prefix}_correlation_matrix.csv")
    print(f"  3. {args.output_prefix}_dendrogram.png")
    print(f"  4. {args.output_prefix}_correlation_heatmap.png")
    
    # Mention pairwise comparison files if they were generated
    if len(gene_spectra_dict) <= 4:
        n_pairwise = len(gene_spectra_dict) * (len(gene_spectra_dict) - 1) // 2
        print(f"  5. {args.output_prefix}_pairwise_k*_vs_k*.png ({n_pairwise} files)")
        print(f"  6. {args.output_prefix}_similar_geps_threshold*.csv (multiple thresholds)")
        print(f"  7. {args.output_prefix}_analysis_summary.txt")
    else:
        print(f"  5. {args.output_prefix}_similar_geps_threshold*.csv (multiple thresholds)")
        print(f"  6. {args.output_prefix}_analysis_summary.txt")
    
    print("\nSee summary file for interpretation guidance!")
    print("="*80)


if __name__ == "__main__":
    main()#!/usr/bin/env python3
"""
cNMF K-Value Comparison: Hierarchical Clustering Analysis
==========================================================

This script analyzes gene expression programs (GEPs) identified by cNMF across
different k values for a SINGLE dataset to identify:
- Which GEPs are consistently found across different k values
- How GEPs cluster together across k values
- Correlation structure of GEPs across k values

Author: Analysis pipeline for cNMF consensus comparison
"""

import argparse
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from scipy.stats import zscore
import warnings
warnings.filterwarnings('ignore')

# Try to import cnmf
try:
    import cnmf
except ImportError:
    print("ERROR: cnmf package not found. Please install it:")
    print("  pip install cnmf")
    sys.exit(1)


def run_consensus_for_k(cnmf_dir, k_value, density_threshold=0.01):
    """
    Run cNMF consensus for a specific k value.
    
    Parameters:
    -----------
    cnmf_dir : str or Path
        Path to the cNMF results directory
    k_value : int
        Number of components (k)
    density_threshold : float
        Density threshold for consensus (default: 0.01)
    
    Returns:
    --------
    str : Path to the gene spectra score file
    """
    print(f"  Running consensus for k={k_value}...")
    
    # Initialize cNMF object
    cnmf_obj = cnmf.cNMF(output_dir=str(cnmf_dir), name=Path(cnmf_dir).name)
    
    # Run consensus
    try:
        cnmf_obj.consensus(k=k_value, 
                          density_threshold=density_threshold,
                          show_clustering=False)
        
        # Construct expected output file path
        # NOTE: cNMF replaces decimal point with underscore in filenames
        # e.g., dt_0.01 becomes dt_0_01
        dt_str = f"{density_threshold:.2f}".replace('.', '_')
        gene_spectra_file = os.path.join(
            cnmf_dir,
            f"{Path(cnmf_dir).name}.gene_spectra_score.k_{k_value}.dt_{dt_str}.txt"
        )
        
        if os.path.exists(gene_spectra_file):
            print(f"    ✓ Consensus completed. Output: {os.path.basename(gene_spectra_file)}")
            return gene_spectra_file
        else:
            print(f"    ✗ Expected output file not found: {gene_spectra_file}")
            return None
            
    except Exception as e:
        print(f"    ✗ Error running consensus for k={k_value}: {str(e)}")
        return None


def load_gene_spectra(file_path, k_value, dataset_name):
    """
    Load gene spectra score matrix and add metadata.
    
    Parameters:
    -----------
    file_path : str
        Path to gene spectra score file
    k_value : int
        K value for this matrix
    dataset_name : str
        Name of the dataset
    
    Returns:
    --------
    pd.DataFrame : Gene spectra matrix with renamed indices
    """
    print(f"  Loading gene spectra for k={k_value}...")
    
    try:
        # Load the matrix (rows = GEPs, columns = genes)
        df = pd.read_csv(file_path, sep='\t', index_col=0)
        
        # Rename indices to include k value: k10_GEP1, k10_GEP2, etc.
        df.index = [f"k{k_value}_GEP{i+1}" for i in range(len(df))]
        
        print(f"    ✓ Loaded {len(df)} GEPs × {len(df.columns)} genes")
        return df
        
    except Exception as e:
        print(f"    ✗ Error loading file: {str(e)}")
        return None


def stack_matrices(gene_spectra_dict):
    """
    Stack all gene spectra matrices vertically.
    
    Parameters:
    -----------
    gene_spectra_dict : dict
        Dictionary mapping k values to gene spectra DataFrames
    
    Returns:
    --------
    pd.DataFrame : Stacked matrix with all GEPs from all k values
    """
    print("\nStacking gene spectra matrices...")
    
    # Get all k values in sorted order
    k_values = sorted(gene_spectra_dict.keys())
    
    # Stack matrices vertically
    stacked = pd.concat([gene_spectra_dict[k] for k in k_values], axis=0)
    
    print(f"  ✓ Stacked matrix: {len(stacked)} total GEPs × {len(stacked.columns)} genes")
    print(f"  ✓ K values included: {', '.join(map(str, k_values))}")
    
    # Print breakdown by k
    for k in k_values:
        n_geps = len(gene_spectra_dict[k])
        print(f"    - k={k}: {n_geps} GEPs")
    
    return stacked


def calculate_correlation_matrix(stacked_matrix):
    """
    Calculate Pearson correlation matrix between all GEPs.
    
    Parameters:
    -----------
    stacked_matrix : pd.DataFrame
        Stacked gene spectra matrix
    
    Returns:
    --------
    pd.DataFrame : Pearson correlation matrix (GEP × GEP)
    """
    print("\nCalculating Pearson correlation matrix...")
    
    # Calculate correlation between rows (GEPs)
    # Each row is a GEP, each column is a gene
    corr_matrix = stacked_matrix.T.corr(method='pearson')
    
    print(f"  ✓ Correlation matrix: {corr_matrix.shape[0]} × {corr_matrix.shape[1]}")
    
    # Print some statistics
    # Exclude diagonal (self-correlation = 1.0)
    off_diagonal = corr_matrix.values[~np.eye(corr_matrix.shape[0], dtype=bool)]
    print(f"  ✓ Correlation statistics (excluding diagonal):")
    print(f"    - Mean: {np.mean(off_diagonal):.3f}")
    print(f"    - Median: {np.median(off_diagonal):.3f}")
    print(f"    - Min: {np.min(off_diagonal):.3f}")
    print(f"    - Max: {np.max(off_diagonal):.3f}")
    
    return corr_matrix


def perform_hierarchical_clustering(correlation_matrix, linkage_method='average'):
    """
    Perform hierarchical clustering on correlation matrix.
    
    Parameters:
    -----------
    correlation_matrix : pd.DataFrame
        Pearson correlation matrix
    linkage_method : str
        Linkage method for clustering (default: 'average')
    
    Returns:
    --------
    np.ndarray : Linkage matrix from hierarchical clustering
    """
    print(f"\nPerforming hierarchical clustering (method: {linkage_method})...")
    
    # Convert correlation to distance: distance = 1 - correlation
    # Higher correlation = lower distance
    distance_matrix = 1 - correlation_matrix.values
    
    # Ensure distance matrix is valid (symmetric, non-negative, zero diagonal)
    distance_matrix = np.clip(distance_matrix, 0, 2)  # Correlation ranges from -1 to 1
    np.fill_diagonal(distance_matrix, 0)
    
    # Convert to condensed distance matrix for scipy
    condensed_dist = squareform(distance_matrix, checks=False)
    
    # Perform hierarchical clustering
    linkage_matrix = linkage(condensed_dist, method=linkage_method)
    
    print(f"  ✓ Hierarchical clustering completed")
    
    return linkage_matrix


def plot_dendrogram(linkage_matrix, labels, output_prefix, figsize=(14, 8)):
    """
    Create and save dendrogram visualization.
    
    Parameters:
    -----------
    linkage_matrix : np.ndarray
        Linkage matrix from hierarchical clustering
    labels : list
        Labels for each GEP (e.g., ['k10_GEP1', 'k10_GEP2', ...])
    output_prefix : str
        Prefix for output file
    figsize : tuple
        Figure size (width, height)
    """
    print("\nCreating dendrogram...")
    
    # Extract k values from labels for coloring
    k_values = [int(label.split('_')[0][1:]) for label in labels]
    unique_k = sorted(set(k_values))
    
    # Create color map
    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_k)))
    k_to_color = {k: colors[i] for i, k in enumerate(unique_k)}
    label_colors = [k_to_color[k] for k in k_values]
    
    # Create figure
    plt.figure(figsize=figsize)
    
    # Plot dendrogram
    dendro = dendrogram(
        linkage_matrix,
        labels=labels,
        leaf_rotation=90,
        leaf_font_size=8,
        color_threshold=0
    )
    
    # Color the labels by k value
    ax = plt.gca()
    xlbls = ax.get_xmajorticklabels()
    for lbl in xlbls:
        label_text = lbl.get_text()
        if label_text in labels:
            idx = labels.index(label_text)
            lbl.set_color(label_colors[idx])
    
    plt.title('Hierarchical Clustering of GEPs Across K Values', fontsize=14, fontweight='bold')
    plt.xlabel('GEP ID', fontsize=12)
    plt.ylabel('Distance (1 - Pearson Correlation)', fontsize=12)
    
    # Add legend
    legend_elements = [plt.Line2D([0], [0], color=k_to_color[k], lw=4, label=f'k={k}')
                      for k in unique_k]
    plt.legend(handles=legend_elements, loc='upper right', title='K Value')
    
    plt.tight_layout()
    
    output_file = f"{output_prefix}_dendrogram.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Dendrogram saved: {output_file}")
    plt.close()


def plot_correlation_heatmap(correlation_matrix, output_prefix, figsize=(12, 10)):
    """
    Create and save correlation heatmap with hierarchical clustering.
    Now includes color-coded axes showing k-value for each GEP.
    
    Parameters:
    -----------
    correlation_matrix : pd.DataFrame
        Pearson correlation matrix
    output_prefix : str
        Prefix for output file
    figsize : tuple
        Figure size (width, height)
    """
    print("\nCreating correlation heatmap...")
    
    # Extract k values from row/column labels
    labels = correlation_matrix.index.tolist()
    k_values = [int(label.split('_')[0][1:]) for label in labels]
    unique_k = sorted(set(k_values))
    
    # Create color map for k values
    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_k)))
    k_to_color = {k: colors[i] for i, k in enumerate(unique_k)}
    
    # Create row and column colors for the heatmap
    row_colors = [k_to_color[k] for k in k_values]
    col_colors = [k_to_color[k] for k in k_values]
    
    # Use seaborn clustermap for automatic clustering with color bars
    g = sns.clustermap(
        correlation_matrix,
        cmap='RdBu_r',
        center=0,
        vmin=-1, vmax=1,
        figsize=figsize,
        dendrogram_ratio=0.15,
        cbar_pos=(0.02, 0.83, 0.03, 0.15),
        linewidths=0,
        xticklabels=True,
        yticklabels=True,
        cbar_kws={'label': 'Pearson Correlation'},
        row_colors=row_colors,
        col_colors=col_colors
    )
    
    # Rotate labels
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), rotation=90, fontsize=6)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), rotation=0, fontsize=6)
    
    # Add title
    plt.suptitle('Pearson Correlation Between GEPs Across K Values', 
                 fontsize=14, fontweight='bold', y=0.98)
    
    # Add legend for k-value colors
    # Create legend patches
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=k_to_color[k], label=f'k={k}')
                      for k in unique_k]
    
    # Add legend to the right side of the plot
    g.ax_heatmap.legend(handles=legend_elements, 
                       title='K Value', 
                       bbox_to_anchor=(1.35, 1.0),
                       loc='upper left',
                       frameon=True)
    
    output_file = f"{output_prefix}_correlation_heatmap.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Correlation heatmap saved: {output_file}")
    plt.close()


def plot_pairwise_k_comparison(correlation_matrix, k1, k2, output_prefix):
    """
    Create pairwise comparison heatmap between GEPs from two different k values.
    
    Parameters:
    -----------
    correlation_matrix : pd.DataFrame
        Full Pearson correlation matrix
    k1 : int
        First k value
    k2 : int
        Second k value
    output_prefix : str
        Prefix for output file
    """
    print(f"  Creating pairwise comparison heatmap: k={k1} vs k={k2}...")
    
    # Get GEPs for each k value
    k1_geps = [gep for gep in correlation_matrix.index if gep.startswith(f'k{k1}_')]
    k2_geps = [gep for gep in correlation_matrix.index if gep.startswith(f'k{k2}_')]
    
    if not k1_geps or not k2_geps:
        print(f"    ✗ No GEPs found for k={k1} or k={k2}")
        return
    
    # Extract the submatrix for this pair
    pairwise_corr = correlation_matrix.loc[k1_geps, k2_geps]
    
    # Calculate figure size based on number of GEPs
    n_k1 = len(k1_geps)
    n_k2 = len(k2_geps)
    figsize = (max(8, n_k2 * 0.5), max(6, n_k1 * 0.5))
    
    # Create heatmap
    plt.figure(figsize=figsize)
    
    # Plot heatmap with better visibility
    sns.heatmap(
        pairwise_corr,
        cmap='RdBu_r',
        center=0,
        vmin=-1, vmax=1,
        annot=True if (n_k1 <= 12 and n_k2 <= 12) else False,  # Only annotate if small
        fmt='.2f',
        cbar_kws={'label': 'Pearson Correlation'},
        linewidths=0.5,
        linecolor='gray',
        square=False
    )
    
    plt.title(f'Pairwise GEP Correlation: k={k1} vs k={k2}', 
             fontsize=14, fontweight='bold', pad=20)
    plt.xlabel(f'GEPs from k={k2}', fontsize=12, fontweight='bold')
    plt.ylabel(f'GEPs from k={k1}', fontsize=12, fontweight='bold')
    plt.tight_layout()
    
    output_file = f"{output_prefix}_pairwise_k{k1}_vs_k{k2}.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"    ✓ Pairwise heatmap saved: {output_file}")
    plt.close()


def create_all_pairwise_comparisons(correlation_matrix, k_values, output_prefix):
    """
    Create pairwise comparison heatmaps for all combinations of k values.
    Only called when there are 4 or fewer k values.
    
    Parameters:
    -----------
    correlation_matrix : pd.DataFrame
        Full Pearson correlation matrix
    k_values : list
        List of k values to compare
    output_prefix : str
        Prefix for output files
    """
    print("\nCreating pairwise k-value comparison heatmaps...")
    print(f"  Number of k values: {len(k_values)}")
    
    k_sorted = sorted(k_values)
    n_comparisons = 0
    
    # Generate all pairwise combinations
    for i in range(len(k_sorted)):
        for j in range(i + 1, len(k_sorted)):
            k1 = k_sorted[i]
            k2 = k_sorted[j]
            plot_pairwise_k_comparison(correlation_matrix, k1, k2, output_prefix)
            n_comparisons += 1
    
    print(f"\n  ✓ Created {n_comparisons} pairwise comparison heatmaps")
    print(f"  ✓ Files: {output_prefix}_pairwise_k*_vs_k*.png")



def identify_similar_geps(correlation_matrix, threshold=0.8, output_prefix=None):
    """
    Identify highly correlated GEPs across different k values.
    
    Parameters:
    -----------
    correlation_matrix : pd.DataFrame
        Pearson correlation matrix
    threshold : float
        Correlation threshold for identifying similar GEPs (default: 0.8)
    output_prefix : str, optional
        If provided, save results to file
    
    Returns:
    --------
    pd.DataFrame : Table of highly correlated GEP pairs
    """
    print(f"\nIdentifying similar GEPs (correlation ≥ {threshold})...")
    
    similar_pairs = []
    
    # Get upper triangle indices (avoid duplicates and self-correlations)
    rows, cols = np.triu_indices_from(correlation_matrix, k=1)
    
    for i, j in zip(rows, cols):
        gep1 = correlation_matrix.index[i]
        gep2 = correlation_matrix.index[j]
        corr = correlation_matrix.iloc[i, j]
        
        if corr >= threshold:
            # Extract k values
            k1 = int(gep1.split('_')[0][1:])
            k2 = int(gep2.split('_')[0][1:])
            
            similar_pairs.append({
                'GEP_1': gep1,
                'GEP_2': gep2,
                'K_1': k1,
                'K_2': k2,
                'Correlation': corr,
                'Cross_K': k1 != k2
            })
    
    if similar_pairs:
        df = pd.DataFrame(similar_pairs)
        df = df.sort_values('Correlation', ascending=False)
        
        print(f"  ✓ Found {len(df)} GEP pairs with correlation ≥ {threshold}")
        
        # Count cross-k vs same-k pairs
        cross_k = df['Cross_K'].sum()
        same_k = len(df) - cross_k
        print(f"    - Cross-k pairs (GEPs from different k): {cross_k}")
        print(f"    - Same-k pairs (GEPs from same k): {same_k}")
        
        if output_prefix:
            output_file = f"{output_prefix}_similar_geps_threshold{threshold:.2f}.csv"
            df.to_csv(output_file, index=False)
            print(f"  ✓ Similar GEPs table saved: {output_file}")
        
        return df
    else:
        print(f"  ✗ No GEP pairs found with correlation ≥ {threshold}")
        return pd.DataFrame()


def create_summary_report(stacked_matrix, correlation_matrix, gene_spectra_dict, 
                         linkage_matrix, output_prefix, args):
    """
    Create a comprehensive text summary of the analysis.
    
    Parameters:
    -----------
    stacked_matrix : pd.DataFrame
        Stacked gene spectra matrix
    correlation_matrix : pd.DataFrame
        Pearson correlation matrix
    gene_spectra_dict : dict
        Dictionary of gene spectra matrices by k value
    linkage_matrix : np.ndarray
        Hierarchical clustering linkage matrix
    output_prefix : str
        Prefix for output file
    args : argparse.Namespace
        Command line arguments
    """
    print("\nGenerating summary report...")
    
    output_file = f"{output_prefix}_analysis_summary.txt"
    
    with open(output_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("cNMF K-VALUE COMPARISON: HIERARCHICAL CLUSTERING ANALYSIS\n")
        f.write("="*80 + "\n\n")
        
        # Analysis parameters
        f.write("ANALYSIS PARAMETERS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Dataset directory: {args.cnmf_dir}\n")
        f.write(f"Dataset name: {Path(args.cnmf_dir).name}\n")
        f.write(f"K values analyzed: {', '.join(map(str, sorted(gene_spectra_dict.keys())))}\n")
        f.write(f"Density threshold: {args.density_threshold}\n")
        f.write(f"Linkage method: {args.linkage_method}\n")
        f.write(f"Number of genes: {stacked_matrix.shape[1]}\n\n")
        
        # GEP counts
        f.write("GEP COUNTS BY K VALUE\n")
        f.write("-" * 80 + "\n")
        total_geps = 0
        for k in sorted(gene_spectra_dict.keys()):
            n_geps = len(gene_spectra_dict[k])
            total_geps += n_geps
            f.write(f"  k={k:3d}: {n_geps:3d} GEPs\n")
        f.write(f"  Total: {total_geps:3d} GEPs\n\n")
        
        # Correlation statistics
        f.write("CORRELATION STATISTICS\n")
        f.write("-" * 80 + "\n")
        off_diag = correlation_matrix.values[~np.eye(correlation_matrix.shape[0], dtype=bool)]
        f.write(f"  Mean correlation: {np.mean(off_diag):.4f}\n")
        f.write(f"  Median correlation: {np.median(off_diag):.4f}\n")
        f.write(f"  Std deviation: {np.std(off_diag):.4f}\n")
        f.write(f"  Min correlation: {np.min(off_diag):.4f}\n")
        f.write(f"  Max correlation: {np.max(off_diag):.4f}\n\n")
        
        # Correlation distribution
        f.write("CORRELATION DISTRIBUTION\n")
        f.write("-" * 80 + "\n")
        bins = [(-1.0, -0.5), (-0.5, 0.0), (0.0, 0.3), (0.3, 0.5), 
                (0.5, 0.7), (0.7, 0.9), (0.9, 1.0)]
        for low, high in bins:
            count = np.sum((off_diag >= low) & (off_diag < high))
            pct = 100 * count / len(off_diag)
            f.write(f"  [{low:5.2f}, {high:5.2f}): {count:6d} pairs ({pct:5.2f}%)\n")
        f.write("\n")
        
        # High correlation pairs analysis
        f.write("HIGHLY CORRELATED GEP PAIRS (correlation ≥ 0.8)\n")
        f.write("-" * 80 + "\n")
        
        similar_df = identify_similar_geps(correlation_matrix, threshold=0.8, 
                                          output_prefix=None)
        
        if len(similar_df) > 0:
            cross_k = similar_df['Cross_K'].sum()
            same_k = len(similar_df) - cross_k
            
            f.write(f"  Total highly correlated pairs: {len(similar_df)}\n")
            f.write(f"  Cross-k pairs (different k values): {cross_k}\n")
            f.write(f"  Same-k pairs (same k value): {same_k}\n\n")
            
            if cross_k > 0:
                f.write("  Cross-k highly correlated pairs (top 20):\n")
                cross_k_pairs = similar_df[similar_df['Cross_K']].head(20)
                for idx, row in cross_k_pairs.iterrows():
                    f.write(f"    {row['GEP_1']:15s} ↔ {row['GEP_2']:15s}  "
                           f"(r = {row['Correlation']:.4f})\n")
        else:
            f.write("  No highly correlated pairs found (threshold = 0.8)\n")
        f.write("\n")
        
        # Interpretation guide
        f.write("="*80 + "\n")
        f.write("INTERPRETATION GUIDE\n")
        f.write("="*80 + "\n\n")
        
        f.write("WHAT DO THESE RESULTS MEAN?\n")
        f.write("-" * 80 + "\n")
        f.write("This analysis helps you understand:\n\n")
        
        f.write("1. CONSISTENCY ACROSS K VALUES:\n")
        f.write("   • High cross-k correlations indicate GEPs that are consistently\n")
        f.write("     identified regardless of the k parameter choice\n")
        f.write("   • These are likely robust, biologically meaningful programs\n\n")
        
        f.write("2. K VALUE STABILITY:\n")
        f.write("   • If GEPs from adjacent k values (e.g., k=10 and k=11) cluster\n")
        f.write("     together, it suggests those k values capture similar biology\n")
        f.write("   • Large jumps in GEP composition between k values may indicate\n")
        f.write("     a transition in model complexity\n\n")
        
        f.write("3. OPTIMAL K SELECTION:\n")
        f.write("   • Look for k values where:\n")
        f.write("     - GEPs are distinct from each other (low intra-k correlation)\n")
        f.write("     - But robust programs are preserved across k (high cross-k correlation)\n")
        f.write("   • The dendrogram can help visualize when new GEPs emerge vs when\n")
        f.write("     existing GEPs split into sub-programs\n\n")
        
        f.write("HOW TO USE THE OUTPUT FILES:\n")
        f.write("-" * 80 + "\n")
        f.write("• Dendrogram: Shows hierarchical relationships between ALL GEPs\n")
        f.write("  - GEPs are colored by k value\n")
        f.write("  - Mixed-color clusters = programs stable across k values\n")
        f.write("  - Single-color clusters = k-specific programs\n\n")
        
        f.write("• Correlation Heatmap: Shows pairwise correlations\n")
        f.write("  - Dark red blocks = highly correlated GEPs\n")
        f.write("  - Off-diagonal blocks = cross-k correlations\n")
        f.write("  - Use to identify which specific GEPs are similar\n\n")
        
        f.write("• Similar GEPs Table: Lists all high-correlation pairs\n")
        f.write("  - Focus on cross-k pairs to find stable programs\n")
        f.write("  - These GEPs likely represent the same biological process\n\n")
        
        f.write("• Stacked Matrix: Full gene weights for all GEPs\n")
        f.write("  - Use for downstream analysis (e.g., gene set enrichment)\n")
        f.write("  - Compare gene loadings for similar GEPs across k values\n\n")
        
        f.write("="*80 + "\n")
        f.write("NEXT STEPS\n")
        f.write("="*80 + "\n")
        f.write("1. Examine the dendrogram and heatmap to identify stable GEPs\n")
        f.write("2. For highly correlated cross-k GEPs, compare their top genes\n")
        f.write("3. Consider the k value that:\n")
        f.write("   a) Captures the key stable programs (seen across multiple k)\n")
        f.write("   b) Provides useful resolution of biological processes\n")
        f.write("   c) Doesn't over-split programs into redundant components\n")
        f.write("4. Use biological knowledge and gene set enrichment to validate\n")
        f.write("   that the selected k captures meaningful programs\n")
        f.write("="*80 + "\n")
    
    print(f"  ✓ Summary report saved: {output_file}")


def main():
    """Main analysis pipeline."""
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Hierarchical clustering analysis of cNMF GEPs across k values',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python cnmf_k_comparison_clustering.py \\
      --cnmf_dir /path/to/dataset \\
      --k_values 8 9 10 11 12 \\
      --density_threshold 0.01 \\
      --output_prefix dataset_k_comparison

This will:
  1. Run cNMF consensus for k=8,9,10,11,12
  2. Load gene spectra score matrices for each k
  3. Stack all matrices vertically
  4. Calculate Pearson correlation between all GEPs
  5. Perform hierarchical clustering
  6. Generate visualizations and summary report
        """
    )
    
    parser.add_argument('--cnmf_dir', required=True,
                       help='Path to cNMF results directory for the dataset')
    parser.add_argument('--k_values', required=True, nargs='+', type=int,
                       help='K values to analyze (space-separated)')
    parser.add_argument('--density_threshold', type=float, default=0.01,
                       help='Density threshold for cNMF consensus (default: 0.01)')
    parser.add_argument('--linkage_method', default='average',
                       choices=['average', 'complete', 'single', 'ward', 'weighted'],
                       help='Linkage method for hierarchical clustering (default: average)')
    parser.add_argument('--output_prefix', default='k_comparison',
                       help='Prefix for output files (default: k_comparison)')
    parser.add_argument('--skip_consensus', action='store_true',
                       help='Skip running consensus (assumes already completed)')
    
    args = parser.parse_args()
    
    # Validate inputs
    cnmf_dir = Path(args.cnmf_dir)
    if not cnmf_dir.exists():
        print(f"ERROR: Directory not found: {cnmf_dir}")
        sys.exit(1)
    
    if len(args.k_values) == 0:
        print("ERROR: Must provide at least one k value")
        sys.exit(1)
    
    # Print header
    print("="*80)
    print("cNMF K-VALUE COMPARISON: HIERARCHICAL CLUSTERING ANALYSIS")
    print("="*80)
    print(f"\nDataset: {cnmf_dir.name}")
    print(f"K values: {', '.join(map(str, sorted(args.k_values)))}")
    print(f"Density threshold: {args.density_threshold}")
    print(f"Linkage method: {args.linkage_method}")
    print("\n" + "="*80)
    
    # Step 1: Run consensus for each k value (unless skipped)
    gene_spectra_files = {}
    
    if not args.skip_consensus:
        print("\nSTEP 1: Running cNMF consensus for each k value")
        print("-"*80)
        
        for k in sorted(args.k_values):
            gene_spectra_file = run_consensus_for_k(
                cnmf_dir, k, args.density_threshold
            )
            if gene_spectra_file:
                gene_spectra_files[k] = gene_spectra_file
            else:
                print(f"  WARNING: Consensus failed for k={k}, skipping...")
    else:
        print("\nSTEP 1: Skipping consensus (--skip_consensus flag set)")
        print("-"*80)
        print("Looking for existing gene spectra files...")
        
        dataset_name = cnmf_dir.name
        # NOTE: cNMF replaces decimal point with underscore in filenames
        dt_str = f"{args.density_threshold:.2f}".replace('.', '_')
        
        for k in sorted(args.k_values):
            expected_file = cnmf_dir / f"{dataset_name}.gene_spectra_score.k_{k}.dt_{dt_str}.txt"
            if expected_file.exists():
                gene_spectra_files[k] = str(expected_file)
                print(f"  ✓ Found k={k}: {expected_file.name}")
            else:
                print(f"  ✗ Not found k={k}: {expected_file.name}")
    
    if not gene_spectra_files:
        print("\nERROR: No gene spectra files found. Cannot proceed.")
        sys.exit(1)
    
    # Step 2: Load gene spectra matrices
    print("\n" + "="*80)
    print("STEP 2: Loading gene spectra matrices")
    print("-"*80)
    
    gene_spectra_dict = {}
    dataset_name = cnmf_dir.name
    
    for k, file_path in sorted(gene_spectra_files.items()):
        df = load_gene_spectra(file_path, k, dataset_name)
        if df is not None:
            gene_spectra_dict[k] = df
    
    if not gene_spectra_dict:
        print("\nERROR: Failed to load any gene spectra matrices. Cannot proceed.")
        sys.exit(1)
    
    # Step 3: Stack matrices
    print("\n" + "="*80)
    print("STEP 3: Stacking gene spectra matrices")
    print("-"*80)
    
    stacked_matrix = stack_matrices(gene_spectra_dict)
    
    # Save stacked matrix
    output_file = f"{args.output_prefix}_stacked_gep_matrix.csv"
    stacked_matrix.to_csv(output_file)
    print(f"\n  ✓ Stacked matrix saved: {output_file}")
    
    # Step 4: Calculate correlation matrix
    print("\n" + "="*80)
    print("STEP 4: Calculating Pearson correlation matrix")
    print("-"*80)
    
    correlation_matrix = calculate_correlation_matrix(stacked_matrix)
    
    # Save correlation matrix
    output_file = f"{args.output_prefix}_correlation_matrix.csv"
    correlation_matrix.to_csv(output_file)
    print(f"\n  ✓ Correlation matrix saved: {output_file}")
    
    # Step 5: Hierarchical clustering
    print("\n" + "="*80)
    print("STEP 5: Performing hierarchical clustering")
    print("-"*80)
    
    linkage_matrix = perform_hierarchical_clustering(
        correlation_matrix, 
        args.linkage_method
    )
    
    # Step 6: Create visualizations
    print("\n" + "="*80)
    print("STEP 6: Creating visualizations")
    print("-"*80)
    
    plot_dendrogram(
        linkage_matrix,
        list(stacked_matrix.index),
        args.output_prefix
    )
    
    plot_correlation_heatmap(
        correlation_matrix,
        args.output_prefix
    )
    
    # Create pairwise comparison heatmaps if 4 or fewer k values
    if len(gene_spectra_dict) <= 4:
        print("\n" + "-"*80)
        print("Creating pairwise k-value comparison heatmaps...")
        print(f"(Generated because number of k values [{len(gene_spectra_dict)}] ≤ 4)")
        print("-"*80)
        create_all_pairwise_comparisons(
            correlation_matrix,
            list(gene_spectra_dict.keys()),
            args.output_prefix
        )
    else:
        print("\n" + "-"*80)
        print(f"Skipping pairwise comparison heatmaps (number of k values [{len(gene_spectra_dict)}] > 4)")
        print("These are only generated when comparing 4 or fewer k values")
        print("-"*80)
    
    # Step 7: Identify similar GEPs
    print("\n" + "="*80)
    print("STEP 7: Identifying similar GEPs")
    print("-"*80)
    
    for threshold in [0.9, 0.8, 0.7]:
        identify_similar_geps(
            correlation_matrix,
            threshold=threshold,
            output_prefix=args.output_prefix
        )
    
    # Step 8: Create summary report
    print("\n" + "="*80)
    print("STEP 8: Creating summary report")
    print("-"*80)
    
    create_summary_report(
        stacked_matrix,
        correlation_matrix,
        gene_spectra_dict,
        linkage_matrix,
        args.output_prefix,
        args
    )
    
    # Final summary
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80)
    print("\nOutput files generated:")
    print(f"  1. {args.output_prefix}_stacked_gep_matrix.csv")
    print(f"  2. {args.output_prefix}_correlation_matrix.csv")
    print(f"  3. {args.output_prefix}_dendrogram.png")
    print(f"  4. {args.output_prefix}_correlation_heatmap.png")
    
    # Mention pairwise comparison files if they were generated
    if len(gene_spectra_dict) <= 4:
        n_pairwise = len(gene_spectra_dict) * (len(gene_spectra_dict) - 1) // 2
        print(f"  5. {args.output_prefix}_pairwise_k*_vs_k*.png ({n_pairwise} files)")
        print(f"  6. {args.output_prefix}_similar_geps_threshold*.csv (multiple thresholds)")
        print(f"  7. {args.output_prefix}_analysis_summary.txt")
    else:
        print(f"  5. {args.output_prefix}_similar_geps_threshold*.csv (multiple thresholds)")
        print(f"  6. {args.output_prefix}_analysis_summary.txt")
    
    print("\nSee summary file for interpretation guidance!")
    print("="*80)


if __name__ == "__main__":
    main()