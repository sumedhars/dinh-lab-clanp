#!/usr/bin/env python3
import os
import sys
import pandas as pd
import scanpy as sc

def main(h5ad_dir, gene_csv, output_txt):
    # --- collect all h5ad files ---
    print("i am here 2")
    h5ad_files = [os.path.join(h5ad_dir, f) for f in os.listdir(h5ad_dir) if f.endswith(".h5ad")]
    if not h5ad_files:
        print(f"No .h5ad files found in {h5ad_dir}")
        sys.exit(1)

    print(f"Found {len(h5ad_files)} .h5ad files.")
    
    # --- get intersection of genes from all .h5ad files ---
    common_genes = None
    for f in h5ad_files:
        print(f"Reading {f} ...")
        adata = sc.read_h5ad(f)
        genes = set(adata.var_names)
        if common_genes is None:
            common_genes = genes
        else:
            common_genes &= genes
        print(f"Current intersection size: {len(common_genes)}")

    print(f"Final intersection across all .h5ad files: {len(common_genes)} genes")

    # --- read gene list from CSV under column 'x' ---
    csv_df = pd.read_csv(gene_csv)
    if "x" not in csv_df.columns:
        print(f"Error: CSV must contain a column named 'x'")
        sys.exit(1)

    csv_genes = csv_df["x"].dropna().astype(str).unique().tolist()
    print(f"Loaded {len(csv_genes)} unique genes from CSV")

    # --- intersection with CSV list ---
    final_genes = sorted(list(common_genes.intersection(csv_genes)))
    print(f"Final intersection size (h5ad ∩ csv): {len(final_genes)} genes")

    # --- write output ---
    os.makedirs(os.path.dirname(output_txt), exist_ok=True)
    with open(output_txt, "w") as f:
        for g in final_genes:
            f.write(g + "\n")

    print(f"Output written to {output_txt}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python find_common_genes.py <h5ad_folder> <gene_csv> <output_txt>")
        sys.exit(1)

    print("i am here 1")
    h5ad_dir = sys.argv[1]
    gene_csv = sys.argv[2]
    output_txt = sys.argv[3]

    main(h5ad_dir, gene_csv, output_txt)
