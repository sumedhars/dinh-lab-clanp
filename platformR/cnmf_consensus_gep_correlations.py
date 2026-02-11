#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Run cNMF consensus across dataset subfolders and compute **GEP×GEP** Pearson
correlations between dataset pairs (one heatmap per dataset pair).

Changes vs original:
- Correlation is computed **between every GEP in dataset_i and every GEP in dataset_j**.
- Outputs a long-form CSV with one row per (dataset_i, dataset_j, gep_i, gep_j).
- Optional: a **pairwise** heatmap per dataset pair (matrix is GEP_i × GEP_j).
- Robust to DataFrame/NumPy inputs; handles NaNs and constant vectors.

Assumptions:
- Every *immediate* subfolder inside --parent_dir is a single dataset run.
- The dataset/run name is exactly that subfolder's basename.
- Each dataset folder (one level below parent) contains:
    ./cnmf_tmp/
    ./[dataset].k_selection_stats.df.npz
    ./[dataset].k_selection.png
    ./[dataset].overdispersed_genes.txt
  (typo variant .overdispered_genes.txt also tolerated)
- Consensus/usage matrix shape from cNMF.load_results is (cells × programs).
- No recursion is performed; deeper nested directories are ignored.
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
matplotlib.use('Agg')  # non-interactive backend (clusters/servers)
import matplotlib.pyplot as plt
import seaborn as sns

# ------------------------------- config knobs ------------------------------- #
# K_RANGE = np.arange(6, 21) # candidate K values (6..20)
NUM_ITER = 200
NUM_HIGHVAR_GENES = 5000
SEED = 14


# ------------------------------- CLI parsing -------------------------------- #
def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run cNMF consensus and compute dataset-pair GEP×GEP correlations "
                    "using single-level dataset discovery (no recursion)."
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
        "--strict_cell_alignment", action="store_true",
        help="If set, require equal #cells across datasets; "
             "otherwise truncate to the minimum for each pair."
    )
    parser.add_argument(
        "--plot_pair_heatmaps",
        action="store_true",
        help="If set, save a GEP×GEP heatmap for each dataset pair."
    )
    args = parser.parse_args()
    return args


# --------------------------------- logging ---------------------------------- #
def log_message(message: str) -> None:
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")
    sys.stdout.flush()


# ----------------------------- folder discovery ----------------------------- #
def immediate_dataset_dirs(parent_dir: Path, ignore_hidden: bool) -> list[Path]:
    """Return sorted list of *immediate* subdirectories of parent_dir. No recursion."""
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
    """Run name is exactly the dataset folder basename (no heuristics)."""
    run_name = dataset_dir.name
    log_message(f" Run name: {run_name}")
    return run_name


def check_expected_artifacts(dataset_dir: Path, run_name: str) -> list[str]:
    """
    Verify expected files/dirs are present. Returns a list of missing items.
    Never crawls beyond the dataset folder level.
    """
    expected = [
        dataset_dir / "cnmf_tmp",
        dataset_dir / f"{run_name}.k_selection_stats.df.npz",
        dataset_dir / f"{run_name}.k_selection.png",
    ]
    # allow either spelling for the overdispersed file
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

    # cNMF locates artifacts by output_dir and name
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

    # Load results for this K and density
    log_message(" Loading results …")
    try:
        usage, spectra_scores, spectra_tpm, top_genes = cnmf_obj.load_results(
            K=selected_k, density_threshold=density_threshold
        )
        consensus = usage  # usage = consensus/mixture matrix (cells × programs)
        log_message(f" Usage (consensus) shape: {consensus.shape} (#cells × #programs)")
        return {
            "name": dataset_dir.name,
            "dataset_dir": str(dataset_dir),
            "run_name": run_name,
            "selected_k": selected_k,
            "usage": usage,
            "spectra_scores": spectra_scores,
            "spectra_tpm": spectra_tpm,
            "top_genes": top_genes,
            "consensus_matrix": consensus,
        }
    except Exception as e:
        log_message(f"ERROR: load_results failed for {dataset_dir.name}: {str(e)}")
        return None


# ------------------------- GEP×GEP correlations (pair) ---------------------- #
def compute_gep_cross_correlations_long(
    results: list[dict],
    strict_cell_alignment: bool
) -> pd.DataFrame:
    """
    For each dataset pair (i, j), compute Pearson correlations between *every*
    GEP in i and *every* GEP in j (rows are GEPs). Return long-form DataFrame:
    ['dataset_i','dataset_j','gep_i','gep_j','r','p'].

    Notes:
    - usage/consensus is cells×programs; we transpose to GEPs×cells.
    - If #cells differ, we truncate both to min cells (unless strict).
    - NaNs/constant vectors -> skip that (i,j,gi,gj).
    """
    names = [r["name"] for r in results]
    mats = [np.asarray(r["consensus_matrix"]).T for r in results]  # (G_i × C_i)

    rows = []
    log_message("Computing GEP×GEP Pearson correlations for each dataset pair …")
    n = len(results)
    for i in range(n):
        for j in range(i + 1, n):
            Mi, Mj = mats[i], mats[j]

            # Align #cells (columns)
            if Mi.shape[1] != Mj.shape[1]:
                msg = (f" Warning: {names[i]} vs {names[j]} have different #cells: "
                       f"{Mi.shape[1]} vs {Mj.shape[1]}")
                if strict_cell_alignment:
                    log_message(msg + " (strict mode: aborting)")
                    raise ValueError("Cell counts differ under --strict_cell_alignment")
                min_cells = min(Mi.shape[1], Mj.shape[1])
                Mi = Mi[:, :min_cells]
                Mj = Mj[:, :min_cells]

            Gi, Gj = Mi.shape[0], Mj.shape[0]

            # Compute all GEP×GEP correlations
            for gi in range(Gi):
                v1 = Mi[gi, :]
                for gj in range(Gj):
                    v2 = Mj[gj, :]

                    mask = np.isfinite(v1) & np.isfinite(v2)
                    if mask.sum() < 2:
                        continue
                    if np.nanstd(v1[mask]) == 0 or np.nanstd(v2[mask]) == 0:
                        continue

                    r, p = pearsonr(v1[mask], v2[mask])
                    if np.isfinite(r) and np.isfinite(p):
                        rows.append({
                            "dataset_i": names[i],
                            "dataset_j": names[j],
                            "gep_i": gi,
                            "gep_j": gj,
                            "r": float(r),
                            "p": float(p),
                        })

    df = pd.DataFrame(rows, columns=["dataset_i", "dataset_j", "gep_i", "gep_j", "r", "p"])
    return df


def save_gep_cross_correlations_long(df: pd.DataFrame, prefix: str) -> Path:
    out = Path(f"{prefix}_gep_cross_correlations.csv")
    df.to_csv(out, index=False)
    log_message(f"GEP×GEP pairwise correlations saved: {out}")
    return out


def plot_dataset_pair_gep_heatmap(
    df_long: pd.DataFrame,
    dataset_i: str,
    dataset_j: str,
    output_png: Path
) -> None:
    """
    Build a GEP_i × GEP_j matrix for a given dataset pair and plot it.
    """
    sub = df_long[(df_long["dataset_i"] == dataset_i) & (df_long["dataset_j"] == dataset_j)]
    if sub.empty:
        log_message(f"No correlations for pair {dataset_i} vs {dataset_j}; skipping heatmap.")
        return

    # Infer matrix extents from observed indices (robust if some pairs are skipped)
    Gi = (int(sub["gep_i"].max()) + 1) if not sub["gep_i"].empty else 0
    Gj = (int(sub["gep_j"].max()) + 1) if not sub["gep_j"].empty else 0
    if Gi == 0 or Gj == 0:
        log_message(f"Cannot infer GEP sizes for {dataset_i} vs {dataset_j}; skipping.")
        return

    mat = np.full((Gi, Gj), np.nan, dtype=float)
    for _, row in sub.iterrows():
        gi, gj, r = int(row["gep_i"]), int(row["gep_j"]), float(row["r"])
        if 0 <= gi < Gi and 0 <= gj < Gj:
            mat[gi, gj] = r

    plt.figure(figsize=(8, 6))
    sns.heatmap(
        mat, vmin=-1, vmax=1, center=0, cmap="coolwarm",
        square=False, linewidths=0.5, cbar_kws={"label": "Pearson r"}
    )
    plt.title(f"{dataset_i} vs {dataset_j}: GEP×GEP Pearson Correlations", fontsize=14, pad=14)
    plt.xlabel(f"GEP (in {dataset_j})")
    plt.ylabel(f"GEP (in {dataset_i})")
    plt.tight_layout()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close()
    log_message(f"Pair heatmap saved: {output_png}")


# --------------------------------- exports ---------------------------------- #
def save_dataset_summary(results: list[dict], prefix: str) -> Path:
    rows = []
    for r in results:
        rows.append({
            "dataset_name": r["name"],
            "dataset_path": r["dataset_dir"],
            "run_name": r["run_name"],
            "selected_k": r["selected_k"],
            "n_cells": r["consensus_matrix"].shape[0],
            "n_programs": r["consensus_matrix"].shape[1],
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
    log_message("cNMF Consensus & Dataset-Pair GEP×GEP Correlation (single-level scan)")
    log_message("=" * 66)
    log_message(f"Parent directory: {parent_dir}")
    log_message(f"Density threshold: {args.density_threshold}")
    log_message(f"Output prefix: {args.output_prefix}")
    log_message(f"Skip consensus: {args.skip_consensus}")
    log_message(f"Ignore hidden: {args.ignore_hidden}")
    log_message(f"Strict align: {args.strict_cell_alignment}")

    # 1) Discover datasets strictly one level down (no recursion)
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

    # 3) Build dataset configs; validate expected artifacts (optionally strict)
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

    # 5) Save dataset summary
    log_message("\n" + "-" * 66)
    log_message("Saving dataset summary")
    log_message("-" * 66)
    _ = save_dataset_summary(results, args.output_prefix)

    # 6) GEP×GEP correlations per dataset pair (no dataset-level aggregates)
    log_message("\n" + "-" * 66)
    log_message("Computing GEP×GEP correlations (per dataset pair)")
    log_message("-" * 66)
    gep_cross_long = compute_gep_cross_correlations_long(
        results, strict_cell_alignment=args.strict_cell_alignment
    )

    # console preview
    if len(gep_cross_long) == 0:
        log_message("WARNING: No valid GEP×GEP correlations were computed (all pairs filtered).")
    else:
        log_message("\nFirst few GEP×GEP correlations:")
        print(gep_cross_long.head().to_string(index=False))

    # 7) Persist long-form CSV
    _ = save_gep_cross_correlations_long(gep_cross_long, args.output_prefix)

    # 8) Optional: heatmap per dataset pair
    if args.plot_pair_heatmaps and len(gep_cross_long) > 0:
        outdir = Path(f"{args.output_prefix}_pair_heatmaps")
        # unique ordered pairs as in computation (i<j)
        pairs = gep_cross_long[["dataset_i", "dataset_j"]].drop_duplicates().values.tolist()
        for dataset_i, dataset_j in pairs:
            plot_dataset_pair_gep_heatmap(
                gep_cross_long,
                dataset_i, dataset_j,
                outdir / f"{dataset_i}__vs__{dataset_j}.heatmap.png"
            )

    # 9) Summary statistics over all GEP×GEP correlations
    log_message("\n" + "-" * 66)
    log_message("Summary statistics (GEP×GEP correlations)")
    log_message("-" * 66)
    if len(gep_cross_long) > 0:
        vals = gep_cross_long["r"].to_numpy()
        log_message(f"Total (dataset pairs × Gi × Gj): {len(vals)}")
        log_message(f"Mean r: {np.mean(vals):.4f}")
        log_message(f"Median r: {np.median(vals):.4f}")
        log_message(f"Min r: {np.min(vals):.4f}")
        log_message(f"Max r: {np.max(vals):.4f}")
        log_message(f"Std dev r: {np.std(vals):.4f}")
    else:
        log_message("No values to summarize.")

    log_message("\nDone.")


if __name__ == "__main__":
    main()
