import argparse, sys, os
import subprocess as sp
import scanpy as sc
import pandas as pd

"""
Run cNMF on a single .h5ad file, selecting the top N HVGs internally.
Uses GNU parallel for the factorization step.

Example command:
python run_parallel_hvg.py --output-dir cnmf_results \
            --name my_analysis --counts data.h5ad \
            --numgenes 5000 \
            -k 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 \
            --n-iter 200 --total-workers 512 --seed 14
"""

def select_hvgs(counts_file, output_dir, name, numgenes):
    """
    Select top N highly variable genes from the .h5ad file.
    Saves a filtered copy and returns (filtered_file_path, actual_num_genes).
    """
    print("=" * 60)
    print("Selecting highly variable genes...")
    print(f"Input file: {counts_file}")
    print(f"Requested HVGs: {numgenes}")

    adata = sc.read_h5ad(counts_file)
    print(f"Original dataset shape: {adata.shape}")

    # Work on a copy so we don't modify the original object unexpectedly
    adata_proc = adata.copy()

    # Normalize + log-transform for HVG selection only
    sc.pp.normalize_total(adata_proc, target_sum=1e4)
    sc.pp.log1p(adata_proc)

    # Select HVGs (flavor='seurat_v3' needs raw counts, so use 'seurat' on log-normed)
    n_top = min(numgenes, adata_proc.n_vars)
    sc.pp.highly_variable_genes(adata_proc, n_top_genes=n_top)
    hvg_mask = adata_proc.var['highly_variable']
    hvg_genes = list(adata_proc.var_names[hvg_mask])
    print(f"HVGs selected: {len(hvg_genes)}")

    # Filter the ORIGINAL (raw counts) adata to these HVGs
    adata_filtered = adata[:, hvg_genes].copy()
    print(f"Filtered dataset shape: {adata_filtered.shape}")

    # Clean index names to prevent cNMF '_index' column issues
    adata_filtered.var.index.name = None
    adata_filtered.obs.index.name = None
    adata_filtered.var_names = adata_filtered.var_names.astype(str)
    adata_filtered.obs_names = adata_filtered.obs_names.astype(str)

    # Save
    filtered_dir = os.path.join(output_dir, name, 'filtered_data')
    os.makedirs(filtered_dir, exist_ok=True)
    filtered_file = os.path.join(filtered_dir, f"{name}_hvg.h5ad")
    adata_filtered.write_h5ad(filtered_file)
    print(f"Saved filtered dataset to: {filtered_file}")
    print("=" * 60)

    return filtered_file, len(hvg_genes)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--name', type=str, required=True,
                        help='Name for this analysis. Output goes to [output-dir]/[name]/...')
    parser.add_argument('--output-dir', type=str, required=True,
                        help='Output directory.')
    parser.add_argument('-c', '--counts', type=str, required=True,
                        help='Input .h5ad counts file (cell x gene, raw counts).')
    parser.add_argument('-k', '--components', type=int, nargs='+',
                        default=list(range(5, 22)),
                        help='Values of k to test. Default: 5-21')
    parser.add_argument('-n', '--n-iter', type=int, default=200,
                        help='Number of NMF iterations per k.')
    parser.add_argument('--total-workers', type=int, default=1,
                        help='Number of GNU parallel workers.')
    parser.add_argument('--seed', type=int, default=None,
                        help='Master random seed.')
    parser.add_argument('--numgenes', type=int, default=5000,
                        help='Number of highly variable genes to select. Default: 5000')
    parser.add_argument('--tpm', type=str, default=None,
                        help='Optional pre-computed TPM matrix (df.npz or txt).')

    args = parser.parse_args()
    argdict = vars(args)

    # ---- HVG selection ----
    filtered_counts_file, num_hvgs = select_hvgs(
        args.counts, args.output_dir, args.name, args.numgenes
    )
    argdict['counts'] = filtered_counts_file
    argdict['numgenes'] = num_hvgs
    print(f"Using {num_hvgs} HVGs for cNMF analysis")

    # Convert components list to space-separated string for CLI
    argdict['components'] = ' '.join([str(k) for k in argdict['components']])

    # Directory containing cnmf.py (assumed to be alongside this script)
    cnmfdir = os.path.dirname(os.path.abspath(sys.argv[0]))
    python_exec = sys.executable

    # ---- prepare ----
    prepare_opts = ['--{} {}'.format(k.replace('_', '-'), argdict[k])
                    for k in argdict if argdict[k] is not None]
    prepare_cmd = f'{python_exec} {cnmfdir}/cnmf.py prepare ' + ' '.join(prepare_opts)
    print(prepare_cmd)
    sp.call(prepare_cmd, shell=True)

    # ---- factorize (GNU parallel) ----
    workind = ' '.join([str(x) for x in range(argdict['total_workers'])])
    factorize_cmd = (
        f'nohup parallel {python_exec} {cnmfdir}/cnmf.py factorize '
        f'--output-dir {argdict["output_dir"]} --name {argdict["name"]} '
        f'--worker-index {{}} ::: {workind}'
    )
    print(factorize_cmd)
    sp.call(factorize_cmd, shell=True)

    # ---- combine ----
    combine_cmd = (
        f'{python_exec} {cnmfdir}/cnmf.py combine '
        f'--output-dir {argdict["output_dir"]} --name {argdict["name"]}'
    )
    print(combine_cmd)
    sp.call(combine_cmd, shell=True)

    # ---- k_selection_plot ----
    kselect_cmd = (
        f'{python_exec} {cnmfdir}/cnmf.py k_selection_plot '
        f'--output-dir {argdict["output_dir"]} --name {argdict["name"]}'
    )
    print(kselect_cmd)
    sp.call(kselect_cmd, shell=True)

    # ---- cleanup iteration files ----
    clean_cmd = f'rm {argdict["output_dir"]}/{argdict["name"]}/cnmf_tmp/*.iter_*.df.npz'
    print(clean_cmd)
    sp.call(clean_cmd, shell=True)


if __name__ == '__main__':
    main()
