#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Run cNMF consensus across dataset subfolders and compute **GEP×GEP** Pearson
correlations **across genes** (one global matrix across all datasets).

What's new (vs previous version that correlated across cells):
- Correlation is computed **per GEP using its gene weights** (spectra), not usage.
- Build a single global correlation matrix of size (sum_i K_i) × (sum_i K_i).
  Example: 4 datasets × 10 GEPs => 40×40 matrix.
- Each correlation is: Pearson( GEP_i (length = #genes), GEP_j (length = #genes) ).
- Aligns spectra by **gene names** across datasets (intersection, consistent order).
  If gene names are missing everywhere, assumes same column order/length.
- Outputs:
    1) Long-form CSV with one row per (dataset_i, gep_i, dataset_j, gep_j, r, p)
    2) Wide CSV with the full correlation matrix (rows/cols annotated)
    3) A single global heatmap PNG (optional block gridlines by dataset)
- Robust to DataFrame/NumPy inputs; handles NaNs and constant vectors.

Assumptions:
- Every *immediate* subfolder inside --parent_dir is a single dataset run.
- The dataset/run name is exactly that subfolder's basename.
- Each dataset folder (one level below parent) contains expected cNMF artifacts.
- cNMF.load_results returns spectra with gene labels (ideally).
"""

from __future__ import annotations
import argparse
import os
from pathlib import Path
from datetime import datetime
import sys
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from cnmf import cNMF

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# ------------------------------- config knobs ------------------------------- #
NUM_ITER = 200
NUM_HIGHVAR_GENES = 5000
SEED = 14


# ------------------------------- CLI parsing -------------------------------- #
def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run cNMF consensus and compute GEP×GEP correlations across genes."
    )
    parser.add_argument(
        "--parent_dir", type=str, required=True,
        help="Parent directory. Each *immediate* subfolder is a dataset run."
    )
    parser.add_argument(
        "--selected_ks", type=str, default=None,
        help='Comma-separated K per dataset, e.g., "10,12,11". '
             "If omitted, --default_k is used for all datasets."
    )
    parser.add_argument(
        "--default_k", type=int, default=10,
        help="Fallback K if --selected_ks not provided or not enough values (default: 10)."
    )
    parser.add_argument(
        "--density_threshold", type=float, default=0.01,
        help="Consensus density threshold (default: 0.01)."
    )
    parser.add_argument(
        "--output_prefix", type=str, default="cnmf_results",
        help="Prefix for output files (CSV, PNG)."
    )
    parser.add_argument(
        "--skip_consensus", action="store_true",
        help="If set, skip running consensus; only load existing results."
    )
    parser.add_argument(
        "--require_expected_artifacts", action="store_true",
        help="If set, require expected cNMF artifacts to exist for each dataset; "
             "otherwise only warn."
    )
    parser.add_argument(
        "--ignore_hidden", action="store_true",
        help="If set, ignore hidden subfolders (those starting with a dot)."
    )
    parser.add_argument(
        "--plot_global_heatmap",
        action="store_true",
        help="If set, save a single global GEP×GEP heatmap across all datasets."
    )
    parser.add_argument(
        "--draw_dataset_blocks",
        action="store_true",
        help="If set with --plot_global_heatmap, draw grid lines separating datasets."
    )
    parser.add_argument(
        "--prefer_spectra",
        choices=["tpm", "scores"], default="tpm",
        help="Which spectra to use if both available (default: tpm)."
    )
    parser.add_argument(
        "--min_common_genes", type=int, default=500,
        help="Minimum # of common genes required across datasets (default: 500)."
    )
    return parser.parse_args()


# --------------------------------- logging ---------------------------------- #
def log_message(message: str) -> None:
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")
    sys.stdout.flush()


# ----------------------------- folder discovery ----------------------------- #
def immediate_dataset_dirs(parent_dir: Path, ignore_hidden: bool) -> list[Path]:
    if not parent_dir.is_dir():
        raise FileNotFoundError(f"Parent directory not found: {parent_dir}")
    subdirs = []
    for child in sorted(parent_dir.iterdir(), key=lambda p: p.name):
        if child.is_dir():
            if ignore_hidden and child.name.startswith("."):
                continue
            subdirs.append(child)
    return subdirs


def run_name_from_folder(dataset_dir: Path) -> str:
    run_name = dataset_dir.name
    log_message(f" Run name: {run_name}")
    return run_name


def check_expected_artifacts(dataset_dir: Path, run_name: str) -> list[str]:
    expected = [
        dataset_dir / "cnmf_tmp",
        dataset_dir / f"{run_name}.k_selection_stats.df.npz",
        dataset_dir / f"{run_name}.k_selection.png",
    ]
    overdisp_ok = any(p.exists() for p in [
        dataset_dir / f"{run_name}.overdispersed_genes.txt",
        dataset_dir / f"{run_name}.overdispered_genes.txt",
    ])
    missing = [str(p) for p in expected if not p.exists()]
    if not overdisp_ok:
        missing.append(str(dataset_dir / f"{run_name}.overdispersed_genes.txt (or .overdispered_genes.txt)"))
    return missing


# ------------------------------ processing core ----------------------------- #
def process_dataset(
    dataset_dir: Path,
    run_name: str,
    selected_k: int,
    density_threshold: float,
    skip_consensus: bool,
) -> dict | None:
    log_message(f"Processing {dataset_dir.name} …")

    cnmf_obj = cNMF(output_dir=str(dataset_dir.parent), name=run_name)

    if not skip_consensus:
        log_message(f" Running consensus: K={selected_k}, density_threshold={density_threshold}")
        try:
            cnmf_obj.consensus(k=selected_k, density_threshold=density_threshold)
            log_message(" Consensus completed")
        except Exception as e:
            log_message(
                " Note: consensus step raised an exception; continuing to load existing results. "
                f"({str(e)[:200]})"
            )
    else:
        log_message(" Skipping consensus (per flag)")

    log_message(" Loading results …")
    try:
        usage, spectra_scores, spectra_tpm, top_genes = cnmf_obj.load_results(
            K=selected_k, density_threshold=density_threshold
        )
        # Shapes can vary by implementation:
        # - usage: cells × programs
        # - spectra_*: often genes × programs (DataFrame) — we will coerce to programs × genes
        return {
            "name": dataset_dir.name,
            "dataset_dir": str(dataset_dir),
            "run_name": run_name,
            "selected_k": selected_k,
            "usage": usage,
            "spectra_scores": spectra_scores,
            "spectra_tpm": spectra_tpm,
        }
    except Exception as e:
        log_message(f"ERROR: load_results failed for {dataset_dir.name}: {str(e)}")
        return None


def _to_programs_by_genes_matrix(
    spectra_obj, prefer: str, dataset_name: str
) -> tuple[np.ndarray, list[str]]:
    """
    Convert spectra (scores or tpm) to a (programs × genes) float array and a gene-name list.
    Accepts pandas DataFrame (preferred) or numpy array (fallback).
    - If DataFrame and shape is genes × programs => transpose to programs × genes.
    - If DataFrame and shape is programs × genes => leave as is.
    - If NumPy without labels => return array and dummy gene names ['g0', 'g1', ...].
    """
    # Choose which spectra container to use
    df_or_arr = None
    if prefer == "tpm" and spectra_obj.get("tpm") is not None:
        df_or_arr = spectra_obj["tpm"]
    elif prefer == "scores" and spectra_obj.get("scores") is not None:
        df_or_arr = spectra_obj["scores"]
    else:
        # fallbacks
        df_or_arr = spectra_obj.get("tpm") or spectra_obj.get("scores")

    if df_or_arr is None:
        raise ValueError(f"{dataset_name}: No spectra_tpm or spectra_scores available.")

    # If it's a pandas object, preserve labels
    try:
        import pandas as pd  # noqa
        if isinstance(df_or_arr, pd.DataFrame):
            df = df_or_arr.copy()
            # Detect orientation by label type: genes should be the index if shape aligns that way
            # Heuristic: if index looks like genes (strings) and columns == K -> assume genes × programs
            # Safer: compare dimensions: smaller dimension is likely programs (K)
            if df.shape[0] >= df.shape[1]:
                # Likely genes × programs -> transpose to programs × genes
                mat = df.T.values.astype(float)
                genes = df.index.astype(str).tolist()
            else:
                # Likely programs × genes already
                mat = df.values.astype(float)
                genes = df.columns.astype(str).tolist()
            return mat, genes
    except Exception:
        pass

    # Otherwise, treat as NumPy array without names
    arr = np.asarray(df_or_arr, dtype=float)
    # We cannot know which axis is programs vs genes; assume programs × genes
    P, G = arr.shape
    genes = [f"g{j}" for j in range(G)]
    return arr, genes


def extract_programs_by_genes(results: list[dict], prefer_spectra: str) -> list[dict]:
    """
    For each dataset, extract a dict:
      {
        'name': str,
        'P': int,
        'genes': [str] or None,
        'M': np.ndarray (programs × genes)
      }
    """
    out = []
    for r in results:
        spectra_dict = {
            "scores": r.get("spectra_scores"),
            "tpm": r.get("spectra_tpm"),
        }
        M, genes = _to_programs_by_genes_matrix(
            spectra_dict, prefer=prefer_spectra, dataset_name=r["name"]
        )
        out.append({
            "name": r["name"],
            "P": M.shape[0],
            "genes": genes,          # may be dummy if unlabeled
            "M": M.astype(float),    # programs × genes
        })
        log_message(f" {r['name']}: spectra matrix (programs × genes) = {M.shape}")
    return out


def align_by_common_genes(pg_list: list[dict], min_common_genes: int) -> tuple[list[dict], list[str]]:
    """
    Intersect gene sets across datasets (if names exist). If *all* are dummy (g0..),
    fall back to assuming same column order and length.
    Returns updated pg_list with matrices restricted to common genes in same order.
    """
    # If any dataset has real gene names (non-g# pattern), try named alignment
    named_lists = []
    for d in pg_list:
        if d["genes"] and not all(g.startswith("g") and g[1:].isdigit() for g in d["genes"]):
            named_lists.append(set(d["genes"]))

    if len(named_lists) > 0:
        common = set.intersection(*named_lists) if named_lists else set()
        common = list(sorted(common))
        if len(common) < min_common_genes:
            raise ValueError(
                f"Only {len(common)} common genes across datasets; "
                f"require at least {min_common_genes}."
            )
        # Reorder each matrix to the common gene order (skip those without names)
        name_to_idx_maps = []
        for d in pg_list:
            name_to_idx = {g: i for i, g in enumerate(d["genes"])}
            name_to_idx_maps.append(name_to_idx)

        aligned = []
        col_indices = [ [name_to_idx_maps[i][g] for g in common if g in name_to_idx_maps[i]] 
                        for i in range(len(pg_list)) ]
        # Ensure all have full coverage
        for i, d in enumerate(pg_list):
            if len(col_indices[i]) != len(common):
                raise ValueError(f"{d['name']}: missing some common genes after alignment.")
            M_aligned = d["M"][:, col_indices[i]]
            aligned.append({**d, "M": M_aligned, "genes": common})
        log_message(f"Aligned by common gene names: {len(common)} genes.")
        return aligned, common

    # Fallback: no real names anywhere — assume identical column order/length
    lengths = [d["M"].shape[1] for d in pg_list]
    if len(set(lengths)) != 1:
        raise ValueError("Gene names unavailable and gene counts differ across datasets.")
    G = lengths[0]
    if G < min_common_genes:
        log_message(f"Warning: only {G} genes available without names; below threshold {min_common_genes}. Proceeding anyway.")
    for d in pg_list:
        d["genes"] = [f"g{j}" for j in range(G)]
    return pg_list, pg_list[0]["genes"]


def build_global_matrix(pg_list: list[dict]) -> tuple[np.ndarray, pd.DataFrame]:
    """
    Stack all programs (rows) across datasets to form S × G matrix where
    S = sum_i P_i. Also build label table with columns:
      ['row_idx','dataset','gep']
    """
    blocks = []
    labels_rows = []
    row_idx = 0
    for d in pg_list:
        P = d["P"]
        blocks.append(d["M"])  # P × G
        for gep in range(P):
            labels_rows.append({"row_idx": row_idx, "dataset": d["name"], "gep": gep})
            row_idx += 1
    X = np.vstack(blocks)  # (S × G)
    labels_df = pd.DataFrame(labels_rows)
    return X, labels_df


def pairwise_pearson_matrix(X: np.ndarray) -> np.ndarray:
    """
    Robust Pearson correlation between rows of X (S × G).
    Handles NaNs and constant rows.
    """
    S, G = X.shape
    R = np.full((S, S), np.nan, dtype=float)
    # Center/scale each row on finite values only
    Xc = X.copy().astype(float)
    for i in range(S):
        mask = np.isfinite(Xc[i, :])
        if mask.sum() < 2:
            continue
        m = np.nanmean(Xc[i, mask])
        s = np.nanstd(Xc[i, mask], ddof=1)
        if s == 0 or not np.isfinite(s):
            continue
        Xc[i, mask] = (Xc[i, mask] - m) / s
        Xc[i, ~mask] = 0.0  # ignore NaNs in dot-products
    # Using dot product since rows are standardized (over non-NaN entries); we normalize by effective N per pair.
    for i in range(S):
        for j in range(i, S):
            mask = np.isfinite(X[i, :]) & np.isfinite(X[j, :])
            n = mask.sum()
            if n < 2:
                continue
            xi = Xc[i, mask]
            xj = Xc[j, mask]
            # After standardization, Pearson r = dot(xi, xj) / (n - 1)
            r = float(np.dot(xi, xj) / (n - 1))
            r = max(min(r, 1.0), -1.0)
            R[i, j] = r
            R[j, i] = r
    return R


def save_global_matrix_csv(R: np.ndarray, labels: pd.DataFrame, prefix: str) -> Path:
    # Build index/columns like 'dataset|gep'
    names = [f"{row.dataset}|{row.gep}" for row in labels.itertuples(index=False)]
    df = pd.DataFrame(R, index=names, columns=names)
    out = Path(f"{prefix}_gep_gene_correlation_matrix.csv")
    df.to_csv(out)
    log_message(f"Global GEP×GEP (across genes) matrix saved: {out}")
    return out


def save_long_form(R: np.ndarray, labels: pd.DataFrame, prefix: str, X: np.ndarray) -> Path:
    # Compute p-values using scipy.stats.pearsonr on the overlapping genes for each pair.
    rows = []
    S = R.shape[0]
    for i in range(S):
        for j in range(i + 1, S):
            r = R[i, j]
            if not np.isfinite(r):
                continue
            # Use same finite mask and pearsonr to get p (redundant but robust)
            mask = np.isfinite(X[i, :]) & np.isfinite(X[j, :])
            if mask.sum() < 2:
                continue
            rr, p = pearsonr(X[i, mask], X[j, mask])
            if not np.isfinite(rr) or not np.isfinite(p):
                continue
            rows.append({
                "dataset_i": labels.loc[i, "dataset"],
                "gep_i": int(labels.loc[i, "gep"]),
                "dataset_j": labels.loc[j, "dataset"],
                "gep_j": int(labels.loc[j, "gep"]),
                "r": float(rr),
                "p": float(p),
            })
    df = pd.DataFrame(rows, columns=["dataset_i","gep_i","dataset_j","gep_j","r","p"])
    out = Path(f"{prefix}_gep_gene_correlations_long.csv")
    df.to_csv(out, index=False)
    log_message(f"Long-form GEP×GEP (across genes) correlations saved: {out}")
    return out


def plot_global_heatmap(R: np.ndarray, labels: pd.DataFrame, prefix: str, draw_blocks: bool) -> None:
    plt.figure(figsize=(10, 8))
    ax = sns.heatmap(
        R, vmin=-1, vmax=1, center=0, cmap="coolwarm",
        square=True, linewidths=0.2, cbar_kws={"label": "Pearson r"}
    )
    # Tick labels
    tick_labels = [f"{row.dataset}|{row.gep}" for row in labels.itertuples(index=False)]
    ax.set_xticks(np.arange(len(tick_labels)) + 0.5)
    ax.set_yticks(np.arange(len(tick_labels)) + 0.5)
    ax.set_xticklabels(tick_labels, rotation=90, fontsize=7)
    ax.set_yticklabels(tick_labels, rotation=0, fontsize=7)
    plt.title("Global GEP×GEP Pearson Correlations (across genes)", fontsize=14, pad=14)
    plt.tight_layout()

    if draw_blocks:
        # Draw grid lines at dataset boundaries
        # Find block sizes by dataset in labels order
        sizes = labels.groupby("dataset").size().tolist()
        cum = np.cumsum(sizes)
        for k in cum[:-1]:
            ax.axhline(k, color="black", lw=1)
            ax.axvline(k, color="black", lw=1)

    outpng = Path(f"{prefix}_global_gep_gene_heatmap.png")
    outpng.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(outpng, dpi=300, bbox_inches="tight")
    plt.close()
    log_message(f"Global heatmap saved: {outpng}")


# --------------------------------- exports ---------------------------------- #
def save_dataset_summary(results: list[dict], prefix: str) -> Path:
    rows = []
    for r in results:
        # usage may be None if load_results changed, so guard shape
        n_cells, n_programs = (np.asarray(r["usage"]).shape if r.get("usage") is not None else (np.nan, np.nan))
        rows.append({
            "dataset_name": r["name"],
            "dataset_path": r["dataset_dir"],
            "run_name": r["run_name"],
            "selected_k": r["selected_k"],
            "n_cells_usage_matrix": n_cells,
            "n_programs_usage_matrix": n_programs,
        })
    df = pd.DataFrame(rows)
    out = Path(f"{prefix}_dataset_summary.csv")
    df.to_csv(out, index=False)
    log_message(f"Dataset summary saved: {out}")
    return out


# ---------------------------------- main ------------------------------------ #
def main() -> None:
    args = parse_arguments()
    parent_dir = Path(args.parent_dir).resolve()

    log_message("=" * 66)
    log_message("cNMF Consensus & Global GEP×GEP Correlation Across Genes")
    log_message("=" * 66)
    log_message(f"Parent directory: {parent_dir}")
    log_message(f"Density threshold: {args.density_threshold}")
    log_message(f"Output prefix: {args.output_prefix}")
    log_message(f"Skip consensus: {args.skip_consensus}")
    log_message(f"Ignore hidden: {args.ignore_hidden}")
    log_message(f"Prefer spectra: {args.prefer_spectra}")
    log_message(f"Min common genes: {args.min_common_genes}")

    # 1) Discover datasets
    log_message("\nDiscovering dataset folders (immediate children only) …")
    dataset_dirs = immediate_dataset_dirs(parent_dir, ignore_hidden=args.ignore_hidden)
    if len(dataset_dirs) == 0:
        log_message("ERROR: No immediate subfolders found; nothing to process.")
        sys.exit(1)
    for d in dataset_dirs:
        log_message(f" - {d}")

    # 2) Parse K choices
    if args.selected_ks:
        provided = [int(k.strip()) for k in args.selected_ks.split(",") if k.strip()]
        if len(provided) != len(dataset_dirs):
            log_message(
                f"Warning: {len(provided)} K values provided for {len(dataset_dirs)} datasets. "
                f"Padding/truncating with default K={args.default_k}."
            )
        if len(provided) < len(dataset_dirs):
            provided.extend([args.default_k] * (len(dataset_dirs) - len(provided)))
        else:
            provided = provided[:len(dataset_dirs)]
        ks = provided
    else:
        ks = [args.default_k] * len(dataset_dirs)
        log_message(f"Using default K={args.default_k} for all datasets.")

    # 3) Validate artifacts
    log_message("\nValidating expected cNMF artifacts per dataset …")
    configs = []
    for dataset_dir, selected_k in zip(dataset_dirs, ks):
        run_name = run_name_from_folder(dataset_dir)
        missing = check_expected_artifacts(dataset_dir, run_name)
        if missing:
            for miss in missing:
                log_message(f" Note: Missing: {miss}")
            if args.require_expected_artifacts:
                log_message(f"ERROR: Missing required artifacts for {dataset_dir.name} (strict mode).")
                sys.exit(1)
        configs.append({
            "dir": dataset_dir,
            "name": dataset_dir.name,
            "run_name": run_name,
            "selected_k": selected_k,
        })
        log_message(f" Ready: {dataset_dir.name} (K={selected_k})")

    # 4) Process datasets
    log_message("\n" + "-" * 66)
    log_message("Processing datasets")
    log_message("-" * 66)
    results = []
    for cfg in configs:
        res = process_dataset(
            dataset_dir=cfg["dir"],
            run_name=cfg["run_name"],
            selected_k=cfg["selected_k"],
            density_threshold=args.density_threshold,
            skip_consensus=args.skip_consensus,
        )
        if res is not None:
            results.append(res)
            log_message(f" OK: {cfg['name']}")
        else:
            log_message(f" FAIL: {cfg['name']}")

    if len(results) < 2:
        log_message("\nERROR: Need at least two successfully processed datasets for correlation analysis.")
        sys.exit(1)

    # 5) Save dataset summary (usage info, if available)
    log_message("\n" + "-" * 66)
    log_message("Saving dataset summary")
    log_message("-" * 66)
    _ = save_dataset_summary(results, args.output_prefix)

    # 6) Extract programs×genes spectra and align genes
    log_message("\n" + "-" * 66)
    log_message("Preparing spectra (programs × genes) and aligning genes")
    log_message("-" * 66)
    pg_list = extract_programs_by_genes(results, prefer_spectra=args.prefer_spectra)
    pg_list_aligned, common_genes = align_by_common_genes(pg_list, args.min_common_genes)

    # 7) Build global S×G matrix (stack GEPs) and compute pairwise Pearson across genes
    log_message("\n" + "-" * 66)
    log_message("Computing global GEP×GEP correlations across genes")
    log_message("-" * 66)
    X, labels = build_global_matrix(pg_list_aligned)
    R = pairwise_pearson_matrix(X)

    # 8) Persist outputs
    _ = save_global_matrix_csv(R, labels, args.output_prefix)
    _ = save_long_form(R, labels, args.output_prefix, X)

    # 9) Optional: plot a single global heatmap
    if args.plot_global_heatmap:
        plot_global_heatmap(R, labels, args.output_prefix, draw_blocks=args.draw_dataset_blocks)

    # 10) Quick summary
    log_message("\n" + "-" * 66)
    log_message("Summary statistics (global GEP×GEP, across genes)")
    log_message("-" * 66)
    vals = R[np.triu_indices_from(R, k=1)]
    vals = vals[np.isfinite(vals)]
    if vals.size > 0:
        log_message(f"Total pairs: {vals.size}")
        log_message(f"Mean r: {np.mean(vals):.4f}")
        log_message(f"Median r: {np.median(vals):.4f}")
        log_message(f"Min r: {np.min(vals):.4f}")
        log_message(f"Max r: {np.max(vals):.4f}")
        log_message(f"Std dev r: {np.std(vals):.4f}")
    else:
        log_message("No valid correlations computed.")

    log_message("\nDone.")


if __name__ == "__main__":
    main()
