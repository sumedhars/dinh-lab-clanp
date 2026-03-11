#!/usr/bin/env python3
"""
cNMF Consensus Runner for Multiple K Values

This script runs cNMF consensus analysis for multiple k values and organizes
the output files into a single directory for easy access.

Usage:
    python cnmf_consensus_runner.py --dataset_dir <path/to/dataset> \
           --k_values <k1> <k2> <k3> --density_threshold <float>
"""

import argparse
import os
import sys
import shutil
from pathlib import Path
from cnmf import cNMF


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Run cNMF consensus for multiple k values and organize outputs",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '--dataset_dir',
        type=str,
        required=True,
        help='Full path to cNMF dataset directory (e.g., /path/to/results/patient1)'
    )
    
    parser.add_argument(
        '--k_values',
        type=int,
        nargs='+',
        required=True,
        help='List of k values to run consensus on (space-separated)'
    )
    
    parser.add_argument(
        '--density_threshold',
        type=float,
        default=0.01,
        help='Density threshold for consensus (default: 0.01)'
    )
    
    parser.add_argument(
        '--output_folder',
        type=str,
        default=None,
        help='Name of folder to create for organized outputs (default: <dataset_name>_consensus_outputs)'
    )
    
    parser.add_argument(
        '--show_clustering',
        action='store_true',
        help='Generate clustering diagnostic plots'
    )
    
    return parser.parse_args()


def verify_cnmf_setup(dataset_dir, k_values):
    """Verify that cNMF has been run through combine step for all k values."""
    print("\n" + "="*70)
    print("VERIFICATION: Checking cNMF setup")
    print("="*70)
    
    dataset_path = Path(dataset_dir)
    name = dataset_path.name
    
    # Check if dataset directory exists
    if not dataset_path.exists():
        print(f"❌ ERROR: Dataset directory not found: {dataset_path}")
        print(f"   Please provide the full path to your cNMF dataset directory")
        sys.exit(1)
    
    print(f"✓ Found dataset directory: {dataset_path}")
    print(f"  Dataset name: {name}")
    
    # Check for cnmf_tmp directory with iteration outputs
    cnmf_tmp = dataset_path / "cnmf_tmp"
    if not cnmf_tmp.exists():
        print(f"⚠ Warning: cnmf_tmp directory not found: {cnmf_tmp}")
        print(f"   This is expected if you're rerunning consensus")
    else:
        print(f"✓ Found iteration outputs directory: {cnmf_tmp}")
    
    # Check for combined results for each k
    missing_k = []
    for k in k_values:
        spectra_pattern = f"{name}.spectra.k_{k}.merged.df.npz"

        # First check in cnmf_tmp/
        spectra_file_tmp = dataset_path / "cnmf_tmp" / spectra_pattern

        # Then fall back to dataset root
        spectra_file_root = dataset_path / spectra_pattern

        if spectra_file_tmp.exists():
            spectra_file = spectra_file_tmp
        elif spectra_file_root.exists():
            spectra_file = spectra_file_root
        else:
            missing_k.append(k)
            continue  # Skip to next k

        print(f"✓ Found merged spectra for k={k}: {spectra_file}")

    
    if missing_k:
        print(f"\n❌ ERROR: Missing combined results for k values: {missing_k}")
        print(f"   Have you run 'cnmf combine' for these k values?")
        print(f"   Expected files like: {name}.spectra.k_<k>.merged.df.npz")
        print(f"   Looking in: {dataset_path}")
        sys.exit(1)
    
    print(f"✓ Found combined results for all k values: {k_values}")
    print("="*70 + "\n")
    
    return True


def run_consensus(dataset_dir, k_values, density_threshold, show_clustering):
    """Run consensus for each k value."""
    print("\n" + "="*70)
    print("RUNNING CONSENSUS")
    print("="*70)
    
    dataset_path = Path(dataset_dir)
    parent_dir = dataset_path.parent
    name = dataset_path.name
    
    print(f"Parent directory (for cNMF): {parent_dir}")
    print(f"Dataset name: {name}")
    
    # Initialize cNMF object with parent directory
    cnmf_obj = cNMF(output_dir=str(parent_dir), name=name)
    
    successful_k = []
    failed_k = []
    
    for k in k_values:
        print(f"\n--- Processing k={k} ---")
        try:
            print(f"Running consensus for k={k} with density_threshold={density_threshold}...")
            cnmf_obj.consensus(
                k=k,
                density_threshold=density_threshold,
                show_clustering=show_clustering,
                refit_usage=False
            )
            print(f"✓ Successfully completed consensus for k={k}")
            successful_k.append(k)
            
        except Exception as e:
            print(f"❌ ERROR: Failed to run consensus for k={k}")
            print(f"   Error message: {str(e)}")
            failed_k.append(k)
    
    print("\n" + "="*70)
    print("CONSENSUS SUMMARY")
    print("="*70)
    print(f"Successfully processed: {successful_k}")
    if failed_k:
        print(f"Failed: {failed_k}")
    print("="*70 + "\n")
    
    return successful_k, failed_k


def organize_outputs(dataset_dir, k_values, density_threshold, output_folder):
    """Create output folder and move consensus files."""
    print("\n" + "="*70)
    print("ORGANIZING OUTPUT FILES")
    print("="*70)
    
    dataset_path = Path(dataset_dir)
    name = dataset_path.name
    
    # Determine output folder name
    if output_folder is None:
        output_folder = f"{name}_consensus_outputs"
    
    # Create output folder in current directory
    output_path = Path.cwd() / output_folder
    output_path.mkdir(exist_ok=True)
    print(f"\n✓ Created/verified output folder: {output_path}")
    
    # Source directory for cNMF consensus outputs (they're in the dataset directory)
    source_dir = dataset_path
    
    # Format density threshold for filenames (0.01 -> 0_01)
    dt_str = str(density_threshold).replace('.', '_')
    
    # File patterns to move
    file_patterns = [
        "{name}.gene_spectra_score.k_{k}.dt_{dt}.txt",
        "{name}.gene_spectra_tpm.k_{k}.dt_{dt}.txt",
        "{name}.usages.k_{k}.dt_{dt}.consensus.txt",
        "{name}.starcat_spectra.k_{k}.dt_{dt}.txt",
        "{name}.spectra.k_{k}.dt_{dt}.consensus.txt"
    ]
    
    moved_files = []
    missing_files = []
    
    for k in k_values:
        print(f"\n--- Moving files for k={k} ---")
        
        for pattern in file_patterns:
            filename = pattern.format(name=name, k=k, dt=dt_str)
            source_file = source_dir / filename
            dest_file = output_path / filename
            
            if source_file.exists():
                shutil.copy2(source_file, dest_file)
                print(f"  ✓ Copied: {filename}")
                moved_files.append(filename)
            else:
                print(f"  ⚠ Missing: {filename}")
                missing_files.append(filename)
        
        # Also copy clustering plots if they exist
        clustering_plot = f"{name}.clustering.k_{k}.dt_{dt_str}.pdf"
        source_plot = source_dir / clustering_plot
        if source_plot.exists():
            shutil.copy2(source_plot, output_path / clustering_plot)
            print(f"  ✓ Copied: {clustering_plot}")
            moved_files.append(clustering_plot)
    
    print("\n" + "="*70)
    print("FILE ORGANIZATION SUMMARY")
    print("="*70)
    print(f"Output folder: {output_path}")
    print(f"Total files copied: {len(moved_files)}")
    if missing_files:
        print(f"Missing files: {len(missing_files)}")
        print("Note: Some files may be missing if consensus failed for certain k values")
    print("="*70 + "\n")
    
    return output_path, moved_files


def create_summary_file(output_path, dataset_name, k_values, density_threshold, moved_files):
    """Create a summary file with run information."""
    summary_file = output_path / "consensus_run_summary.txt"
    
    with open(summary_file, 'w') as f:
        f.write("="*70 + "\n")
        f.write("cNMF CONSENSUS RUN SUMMARY\n")
        f.write("="*70 + "\n\n")
        
        f.write(f"Dataset name: {dataset_name}\n")
        f.write(f"K values processed: {k_values}\n")
        f.write(f"Density threshold: {density_threshold}\n")
        f.write(f"Number of files generated: {len(moved_files)}\n\n")
        
        f.write("="*70 + "\n")
        f.write("OUTPUT FILES\n")
        f.write("="*70 + "\n\n")
        
        # Group files by k value
        for k in k_values:
            f.write(f"K = {k}:\n")
            k_files = [f for f in moved_files if f".k_{k}." in f]
            for fname in sorted(k_files):
                f.write(f"  - {fname}\n")
            f.write("\n")
        
        f.write("="*70 + "\n")
        f.write("FILE DESCRIPTIONS\n")
        f.write("="*70 + "\n\n")
        
        f.write("1. gene_spectra_score files:\n")
        f.write("   Z-score normalized gene expression program (GEP) matrix\n")
        f.write("   Rows = genes, Columns = GEPs\n")
        f.write("   Values = Z-scores indicating gene importance in each program\n\n")
        
        f.write("2. gene_spectra_tpm files:\n")
        f.write("   TPM-normalized gene expression program matrix\n")
        f.write("   Rows = genes, Columns = GEPs\n")
        f.write("   Values = TPM units\n\n")
        
        f.write("3. usages files:\n")
        f.write("   Cell usage matrix (how much each cell uses each GEP)\n")
        f.write("   Rows = cells, Columns = GEPs\n")
        f.write("   Values = normalized usage (sums to 1 for each cell)\n\n")
        
        f.write("4. clustering files (if generated):\n")
        f.write("   Diagnostic clustergram showing GEP similarity\n")
        f.write("   PDF format\n\n")
        
        f.write("="*70 + "\n")
        f.write("NEXT STEPS\n")
        f.write("="*70 + "\n\n")
        
        f.write("1. Load the gene_spectra_score files to examine gene programs\n")
        f.write("2. Load the usages files to analyze program activity across cells\n")
        f.write("3. Perform downstream analysis:\n")
        f.write("   - Gene set enrichment analysis on top genes per GEP\n")
        f.write("   - Compare GEP usage across cell types/conditions\n")
        f.write("   - Correlate GEPs with cell metadata\n")
        f.write("4. To compare GEPs across k values, consider using:\n")
        f.write("   cnmf_k_comparison_clustering.py\n\n")
    
    print(f"✓ Created summary file: {summary_file}")
    return summary_file


def main():
    """Main execution function."""
    args = parse_arguments()
    
    dataset_path = Path(args.dataset_dir)
    dataset_name = dataset_path.name
    
    print("\n" + "="*70)
    print("cNMF CONSENSUS RUNNER")
    print("="*70)
    print(f"Dataset directory: {args.dataset_dir}")
    print(f"Dataset name: {dataset_name}")
    print(f"K values: {args.k_values}")
    print(f"Density threshold: {args.density_threshold}")
    print("="*70)
    
    # Step 1: Verify cNMF setup
    verify_cnmf_setup(args.dataset_dir, args.k_values)
    
    # Step 2: Run consensus for all k values
    successful_k, failed_k = run_consensus(
        args.dataset_dir,
        args.k_values,
        args.density_threshold,
        args.show_clustering
    )
    
    # Step 3: Organize outputs
    if successful_k:
        output_path, moved_files = organize_outputs(
            args.dataset_dir,
            successful_k,
            args.density_threshold,
            args.output_folder
        )
        
        # Step 4: Create summary file
        summary_file = create_summary_file(
            output_path,
            dataset_name,
            successful_k,
            args.density_threshold,
            moved_files
        )
        
        print("\n" + "="*70)
        print("✓ ANALYSIS COMPLETE!")
        print("="*70)
        print(f"\nAll consensus files are in: {output_path}")
        print(f"Read the summary file: {summary_file}")
        print("\nFiles are ready for downstream analysis!")
        print("="*70 + "\n")
        
    else:
        print("\n" + "="*70)
        print("❌ ERROR: No k values were successfully processed")
        print("="*70)
        sys.exit(1)
    
    if failed_k:
        print(f"\n⚠ Warning: Some k values failed: {failed_k}")
        sys.exit(1)


if __name__ == "__main__":
    main()
