import argparse, sys, os
import subprocess as sp
import scanpy as sc
import pandas as pd

"""
Run all of the steps through plotting the K selection plot of cNMF sequentially using GNU
parallel to run the factorization steps in parallel. This version filters the input .h5ad
file to keep only common genes specified in a genes file before running cNMF.

Example command:
python run_parallel_filtered.py --output-dir $output_dir \
            --name test --counts path_to_counts.h5ad \
            --common-genes-file output/common_genes.txt \
            -k 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 --n-iter 100 --total-workers 4 \
            --numgenes 2000 --seed 5
"""

def filter_to_common_genes(counts_file, common_genes_file, output_dir, name):
    """
    Filter the .h5ad file to keep only genes in the common_genes_file.
    Returns a tuple of (filtered_file_path, num_common_genes).
    """
    print("=" * 60)
    print("Filtering dataset to common genes...")
    print(f"Input file: {counts_file}")
    print(f"Common genes file: {common_genes_file}")
    
    # Read the AnnData object
    adata = sc.read_h5ad(counts_file)
    print(f"Original dataset shape: {adata.shape}")
    print(f"Original genes: {adata.n_vars}")
    
    # Read the common genes list
    with open(common_genes_file, 'r') as f:
        common_genes = [line.strip() for line in f if line.strip()]
    print(f"Number of common genes: {len(common_genes)}")
    
    # Filter to common genes that exist in the dataset
    genes_in_dataset = set(adata.var_names)
    common_genes_in_dataset = [g for g in common_genes if g in genes_in_dataset]
    
    if len(common_genes_in_dataset) == 0:
        raise ValueError("No common genes found in the dataset! Check gene name format.")
    
    print(f"Common genes found in dataset: {len(common_genes_in_dataset)}")
    
    # Filter the AnnData object
    adata_filtered = adata[:, common_genes_in_dataset].copy()
    print(f"Filtered dataset shape: {adata_filtered.shape}")
    
    # Create output directory for filtered data
    filtered_dir = os.path.join(output_dir, name, 'filtered_data')
    os.makedirs(filtered_dir, exist_ok=True)
    
    # Save filtered dataset
    filtered_file = os.path.join(filtered_dir, f"{name}_filtered.h5ad")
    adata_filtered.write(filtered_file)
    print(f"Saved filtered dataset to: {filtered_file}")
    print("=" * 60)
    
    return filtered_file, len(common_genes_in_dataset)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--name', type=str, help='[all] Name for this analysis. All output will be placed in [output-dir]/[name]/...', nargs='?', default=None)
    parser.add_argument('--output-dir', type=str, help='[all] Output directory. All output will be placed in [output-dir]/[name]/...', nargs='?')
    parser.add_argument('-c', '--counts', type=str, help='[prepare] Input counts in cell x gene matrix as df.npz, h5ad, or tab separated txt file')
    parser.add_argument('--common-genes-file', type=str, help='[prepare] File containing common genes to filter to, one gene per line', required=True)
    parser.add_argument('-k', '--components', type=int, help='[prepare] Number of components (k) for matrix factorization. Several can be specified with "-k 8 9 10"', nargs='+', default=list(range(5, 22)))  # Default: k=5-21
    parser.add_argument('-n', '--n-iter', type=int, help='[prepare] Number of iterations for each factorization', default=100)
    parser.add_argument('--total-workers', type=int, help='[all] Total workers that are working together.', default=1)
    parser.add_argument('--seed', type=int, help='[prepare] Master seed for generating the seed list.', default=None)
    parser.add_argument('--numgenes', type=int, help='[prepare] Number of high variance genes to use for matrix factorization.', default=2000)  # Default: 2000
    parser.add_argument('--tpm', type=str, help='[prepare] Pre-computed TPM values as df.npz or tab separated txt file. Cell x Gene matrix. If none is provided, TPM will be calculated automatically. This can be helpful if a particular normalization is desired.', default=None)

    # Collect args
    args = parser.parse_args()
    argdict = vars(args)

    # Filter the counts file to common genes
    filtered_counts_file, num_common_genes = filter_to_common_genes(
        args.counts, 
        args.common_genes_file,
        args.output_dir,
        args.name
    )
    
    # Update the counts file to use the filtered version
    argdict['counts'] = filtered_counts_file
    
    # Set numgenes to the number of common genes
    if argdict['numgenes'] != num_common_genes:
        if argdict['numgenes'] is not None:
            print(f"NOTE: Overriding --numgenes {argdict['numgenes']} with {num_common_genes} (number of common genes)")
        argdict['numgenes'] = num_common_genes
    print(f"Using {num_common_genes} genes for cNMF analysis")
    
    # Remove common_genes_file from argdict as it's not a cNMF argument
    del argdict['common_genes_file']

    # convert components from list to string
    argdict['components'] = ' '.join([str(k) for k in argdict['components']])
    
    # Directory containing cNMF and this script
    cnmfdir = os.path.dirname(sys.argv[0])
    if len(cnmfdir) == 0: cnmfdir = '.'   
 
    # Get the Python executable being used to run this script
    python_exec = sys.executable
    
    # Run prepare
    prepare_opts = ['--{} {}'.format(k.replace('_', '-'),argdict[k]) for k in argdict.keys() if argdict[k] is not None]
    prepare_cmd = '{} {}/cnmf.py prepare '.format(python_exec, cnmfdir)
    prepare_cmd += ' '.join(prepare_opts)
    print(prepare_cmd)
    sp.call(prepare_cmd, shell=True)

    # Run factorize
    workind = ' '.join([str(x) for x in range(argdict['total_workers'])])
    factorize_cmd = 'nohup parallel {} {}/cnmf.py factorize --output-dir {} --name {} --worker-index {{}} ::: {}'.format(
        python_exec, cnmfdir, argdict['output_dir'], argdict['name'], workind)
    print(factorize_cmd)
    sp.call(factorize_cmd, shell=True)

    # Run combine
    combine_cmd = '{} {}/cnmf.py combine --output-dir {} --name {}'.format(
        python_exec, cnmfdir, argdict['output_dir'], argdict['name'])
    print(combine_cmd)
    sp.call(combine_cmd, shell=True)

    # Plot K selection
    Kselect_cmd = '{} {}/cnmf.py k_selection_plot --output-dir {} --name {}'.format(
        python_exec, cnmfdir, argdict['output_dir'], argdict['name'])
    print(Kselect_cmd)
    sp.call(Kselect_cmd, shell=True)

    # Delete individual iteration files
    clean_cmd = 'rm %s/%s/cnmf_tmp/*.iter_*.df.npz' % (argdict['output_dir'], argdict['name'])
    print(clean_cmd)
    sp.call(clean_cmd, shell=True)


if __name__ == '__main__':
    main()