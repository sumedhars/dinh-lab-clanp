#!/usr/bin/env python
import argparse, os, sys, glob, shutil
from pathlib import Path
import scanpy as sc

"""
Run cNMF on the union of the top-N highly variable genes computed separately
within each cell type, instead of HVGs computed across the whole dataset.

This is a local/terminal counterpart to run_parallel_hvg.py (which was run on
an HPC cluster via SLURM + GNU parallel). Instead of shelling out to GNU
parallel for the factorize step, it calls cNMF's built-in
factorize_multi_process(), which parallelizes with Python's multiprocessing.Pool
in-process -- no GNU parallel / job scheduler required.

Cell type annotations are read from an .obs column, following the same
pattern used in notebooks/subset_by_annotation.ipynb (e.g. "global.cluster4").

After combine + k_selection_plot, it also runs consensus for each k value and
organizes the outputs into a single folder with a summary file, matching what
cnmf_consensus_runner.py does (pass --skip-consensus to stop before this step).

Example command:
python run_cnmf_celltype_hvg_union.py \
    --output-dir cnmf_results --name my_analysis \
    --counts data/hnc_myeloid_2021.h5ad \
    --celltype-col global.cluster4 \
    --numgenes 1000 \
    -k 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 \
    --n-iter 100 --total-workers 8 --seed 14 \
    --density-threshold 0.01
"""


def select_union_hvgs_per_celltype(counts_file, celltype_col, output_dir, name,
                                    numgenes, min_cells=10):
    """
    For each unique label in adata.obs[celltype_col], compute the top `numgenes`
    highly variable genes within that cell type, then take the union across all
    cell types. Filters the ORIGINAL (raw counts) adata to that union and saves it.

    Returns (filtered_file_path, genes_file_path, num_union_genes).
    """
    print("=" * 60)
    print("Selecting per-cell-type highly variable genes...")
    print(f"Input file: {counts_file}")
    print(f"Cell type column: {celltype_col}")
    print(f"Top HVGs per cell type: {numgenes}")

    adata = sc.read_h5ad(counts_file)
    print(f"Original dataset shape: {adata.shape}")

    if celltype_col not in adata.obs.columns:
        raise ValueError(
            f"'{celltype_col}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    # Normalize + log-transform for HVG selection only (raw counts stay untouched)
    adata_proc = adata.copy()
    sc.pp.normalize_total(adata_proc, target_sum=1e4)
    sc.pp.log1p(adata_proc)

    labels = adata_proc.obs[celltype_col].astype(str)
    cell_types = sorted(labels.unique())
    print(f"Found {len(cell_types)} cell types in '{celltype_col}'\n")

    union_genes = set()
    for ct in cell_types:
        mask = (labels == ct).values
        n_cells = mask.sum()
        if n_cells < min_cells:
            print(f"  {ct}: {n_cells} cells -- skipped (< {min_cells} cells)")
            continue

        sub = adata_proc[mask]
        n_top = min(numgenes, sub.n_vars)
        try:
            hvg_stats = sc.pp.highly_variable_genes(sub, n_top_genes=n_top, inplace=False)
            hvgs = set(sub.var_names[hvg_stats['highly_variable'].values])
        except Exception as e:
            print(f"  {ct}: {n_cells} cells -- HVG selection failed ({e}), skipped")
            continue

        print(f"  {ct}: {n_cells} cells -> {len(hvgs)} HVGs")
        union_genes |= hvgs

    union_genes = sorted(union_genes)
    print(f"\nUnion of HVGs across {len(cell_types)} cell types: {len(union_genes)} genes")

    if len(union_genes) == 0:
        raise RuntimeError("No genes were selected -- check celltype-col and min-cells.")

    # Filter the ORIGINAL (raw counts) adata to the union of genes
    adata_filtered = adata[:, union_genes].copy()
    print(f"Filtered dataset shape: {adata_filtered.shape}")

    # Clean index names to prevent cNMF '_index' column issues
    adata_filtered.var.index.name = None
    adata_filtered.obs.index.name = None
    adata_filtered.var_names = adata_filtered.var_names.astype(str)
    adata_filtered.obs_names = adata_filtered.obs_names.astype(str)

    filtered_dir = os.path.join(output_dir, name, 'filtered_data')
    os.makedirs(filtered_dir, exist_ok=True)

    filtered_file = os.path.join(filtered_dir, f"{name}_hvg_union.h5ad")
    adata_filtered.write_h5ad(filtered_file)
    print(f"Saved filtered dataset to: {filtered_file}")

    genes_file = os.path.join(filtered_dir, f"{name}_hvg_union_genes.txt")
    with open(genes_file, 'w') as f:
        f.write('\n'.join(union_genes))
    print(f"Saved union gene list to: {genes_file}")
    print("=" * 60)

    return filtered_file, genes_file, len(union_genes)


def run_consensus(cnmf_obj, k_values, density_threshold, show_clustering):
    """Run consensus for each k value. Ported from cnmf_consensus_runner.py."""
    print("\n" + "=" * 70)
    print("RUNNING CONSENSUS")
    print("=" * 70)

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
                refit_usage=False,
            )
            print(f"Successfully completed consensus for k={k}")
            successful_k.append(k)
        except Exception as e:
            print(f"ERROR: Failed to run consensus for k={k}: {e}")
            failed_k.append(k)

    print("\n" + "=" * 70)
    print("CONSENSUS SUMMARY")
    print("=" * 70)
    print(f"Successfully processed: {successful_k}")
    if failed_k:
        print(f"Failed: {failed_k}")
    print("=" * 70 + "\n")

    return successful_k, failed_k


def organize_outputs(dataset_dir, name, k_values, density_threshold, output_folder):
    """Create output folder and copy consensus files. Ported from cnmf_consensus_runner.py."""
    print("\n" + "=" * 70)
    print("ORGANIZING OUTPUT FILES")
    print("=" * 70)

    dataset_path = Path(dataset_dir)

    if output_folder is None:
        output_folder = f"{name}_consensus_outputs"

    output_path = Path.cwd() / output_folder
    output_path.mkdir(exist_ok=True)
    print(f"\nCreated/verified output folder: {output_path}")

    dt_str = str(density_threshold).replace('.', '_')

    file_patterns = [
        "{name}.gene_spectra_score.k_{k}.dt_{dt}.txt",
        "{name}.gene_spectra_tpm.k_{k}.dt_{dt}.txt",
        "{name}.usages.k_{k}.dt_{dt}.consensus.txt",
        "{name}.starcat_spectra.k_{k}.dt_{dt}.txt",
        "{name}.spectra.k_{k}.dt_{dt}.consensus.txt",
    ]

    moved_files = []
    missing_files = []

    for k in k_values:
        print(f"\n--- Copying files for k={k} ---")

        for pattern in file_patterns:
            filename = pattern.format(name=name, k=k, dt=dt_str)
            source_file = dataset_path / filename
            dest_file = output_path / filename

            if source_file.exists():
                shutil.copy2(source_file, dest_file)
                print(f"  Copied: {filename}")
                moved_files.append(filename)
            else:
                print(f"  Missing: {filename}")
                missing_files.append(filename)

        clustering_plot = f"{name}.clustering.k_{k}.dt_{dt_str}.png"
        source_plot = dataset_path / clustering_plot
        if source_plot.exists():
            shutil.copy2(source_plot, output_path / clustering_plot)
            print(f"  Copied: {clustering_plot}")
            moved_files.append(clustering_plot)

    print("\n" + "=" * 70)
    print("FILE ORGANIZATION SUMMARY")
    print("=" * 70)
    print(f"Output folder: {output_path}")
    print(f"Total files copied: {len(moved_files)}")
    if missing_files:
        print(f"Missing files: {len(missing_files)}")
        print("Note: Some files may be missing if consensus failed for certain k values")
    print("=" * 70 + "\n")

    return output_path, moved_files


def create_summary_file(output_path, dataset_name, k_values, density_threshold, moved_files):
    """Create a summary file with run information. Ported from cnmf_consensus_runner.py."""
    summary_file = output_path / "consensus_run_summary.txt"

    with open(summary_file, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("cNMF CONSENSUS RUN SUMMARY\n")
        f.write("=" * 70 + "\n\n")

        f.write(f"Dataset name: {dataset_name}\n")
        f.write(f"K values processed: {k_values}\n")
        f.write(f"Density threshold: {density_threshold}\n")
        f.write(f"Number of files generated: {len(moved_files)}\n\n")

        f.write("=" * 70 + "\n")
        f.write("OUTPUT FILES\n")
        f.write("=" * 70 + "\n\n")

        for k in k_values:
            f.write(f"K = {k}:\n")
            k_files = [x for x in moved_files if f".k_{k}." in x]
            for fname in sorted(k_files):
                f.write(f"  - {fname}\n")
            f.write("\n")

        f.write("=" * 70 + "\n")
        f.write("FILE DESCRIPTIONS\n")
        f.write("=" * 70 + "\n\n")

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
        f.write("   Diagnostic clustergram showing GEP similarity\n\n")

        f.write("=" * 70 + "\n")
        f.write("NEXT STEPS\n")
        f.write("=" * 70 + "\n\n")

        f.write("1. Load the gene_spectra_score files to examine gene programs\n")
        f.write("2. Load the usages files to analyze program activity across cells\n")
        f.write("3. Perform downstream analysis:\n")
        f.write("   - Gene set enrichment analysis on top genes per GEP\n")
        f.write("   - Compare GEP usage across cell types/conditions\n")
        f.write("   - Correlate GEPs with cell metadata\n")

    print(f"Created summary file: {summary_file}")
    return summary_file


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--name', type=str, required=True,
                        help='Name for this analysis. Output goes to [output-dir]/[name]/...')
    parser.add_argument('--output-dir', type=str, required=True,
                        help='Output directory.')
    parser.add_argument('-c', '--counts', type=str, required=True,
                        help='Input .h5ad counts file (cell x gene, raw counts).')
    parser.add_argument('--celltype-col', type=str, required=True,
                        help='adata.obs column with cell type / annotation labels '
                             '(e.g. "global.cluster4", see subset_by_annotation.ipynb).')
    parser.add_argument('-k', '--components', type=int, nargs='+',
                        default=list(range(5, 22)),
                        help='Values of k to test. Default: 5-21')
    parser.add_argument('-n', '--n-iter', type=int, default=100,
                        help='Number of NMF iterations per k.')
    parser.add_argument('--total-workers', type=int, default=1,
                        help='Number of local processes for factorization '
                             '(uses Python multiprocessing, no GNU parallel needed).')
    parser.add_argument('--seed', type=int, default=None,
                        help='Master random seed.')
    parser.add_argument('--numgenes', type=int, default=1000,
                        help='Number of top highly variable genes to select PER cell type '
                             'before taking the union. Default: 1000')
    parser.add_argument('--min-cells', type=int, default=10,
                        help='Skip cell types with fewer than this many cells. Default: 10')
    parser.add_argument('--tpm', type=str, default=None,
                        help='Optional pre-computed TPM matrix (df.npz or txt).')
    parser.add_argument('--keep-iterations', action='store_true', default=False,
                        help='Keep per-iteration spectra files after combining (default: delete them).')
    parser.add_argument('--density-threshold', type=float, default=0.01,
                        help='[consensus] Density threshold for consensus filtering. Default: 0.01')
    parser.add_argument('--show-clustering', action='store_true', default=False,
                        help='[consensus] Generate clustering diagnostic plots.')
    parser.add_argument('--consensus-k-values', type=int, nargs='+', default=None,
                        help='[consensus] k values to run consensus for. Defaults to the same '
                             'values passed to -k/--components.')
    parser.add_argument('--consensus-output-folder', type=str, default=None,
                        help='[consensus] Name of folder (created in the current directory) to '
                             'collect organized consensus outputs. Default: <name>_consensus_outputs')
    parser.add_argument('--skip-consensus', action='store_true', default=False,
                        help='Stop after k_selection_plot; skip consensus + output organization.')

    args = parser.parse_args()

    # ---- per-cell-type HVG union selection ----
    filtered_counts_file, genes_file, n_union = select_union_hvgs_per_celltype(
        args.counts, args.celltype_col, args.output_dir, args.name,
        args.numgenes, min_cells=args.min_cells
    )
    print(f"Using {n_union} union HVGs for cNMF analysis\n")

    # Import cnmf.py from alongside this script
    cnmfdir = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, cnmfdir)
    from cnmf import cNMF

    cnmf_obj = cNMF(output_dir=args.output_dir, name=args.name)

    # ---- prepare ----
    # genes_file pins factorization to exactly the union genes we already computed,
    # so cnmf.py won't recompute (and potentially re-narrow) HVGs on its own.
    print("Running prepare...")
    cnmf_obj.prepare(
        counts_fn=filtered_counts_file,
        components=args.components,
        n_iter=args.n_iter,
        tpm_fn=args.tpm,
        seed=args.seed,
        genes_file=genes_file,
        num_highvar_genes=n_union,
    )

    # ---- factorize (local multiprocessing, in-process) ----
    print(f"Running factorization with {args.total_workers} local worker(s)...")
    cnmf_obj.factorize_multi_process(total_workers=args.total_workers)

    # ---- combine ----
    print("Combining iterations...")
    cnmf_obj.combine(components=args.components)

    # ---- k_selection_plot ----
    print("Generating k selection plot...")
    cnmf_obj.k_selection_plot(close_fig=True)

    # ---- cleanup iteration files ----
    if not args.keep_iterations:
        tmp_dir = os.path.join(args.output_dir, args.name, 'cnmf_tmp')
        iter_files = glob.glob(os.path.join(tmp_dir, '*.iter_*.df.npz'))
        print(f"Removing {len(iter_files)} per-iteration spectra file(s)...")
        for f in iter_files:
            os.remove(f)

    # ---- consensus (multiple k values) + output organization ----
    if not args.skip_consensus:
        consensus_k_values = args.consensus_k_values or args.components

        successful_k, failed_k = run_consensus(
            cnmf_obj, consensus_k_values, args.density_threshold, args.show_clustering
        )

        if successful_k:
            dataset_dir = os.path.join(args.output_dir, args.name)
            output_path, moved_files = organize_outputs(
                dataset_dir, args.name, successful_k, args.density_threshold,
                args.consensus_output_folder
            )
            create_summary_file(output_path, args.name, successful_k,
                                 args.density_threshold, moved_files)
            print(f"\nAll consensus files are in: {output_path}")
        else:
            print("WARNING: No k values were successfully processed during consensus")

        if failed_k:
            print(f"WARNING: Some k values failed consensus: {failed_k}")

    print("Done.")


if __name__ == '__main__':
    main()
