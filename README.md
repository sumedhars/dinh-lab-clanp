# CLANP

**Cluster Label Assessment via NMF Programs**

CLANP evaluates scRNA-seq cluster annotations by tying them to consensus NMF (cNMF) gene expression programs across resolutions, then scoring whether each annotated label corresponds to a pure cell type, a mixture that should be split, or two labels that should be merged.

Developed in the [Dinh Lab](https://dinhlab.oncology.wisc.edu/), McArdle Laboratory for Cancer Research, UW–Madison.

---

## What it does

Single-cell clustering is usually done at one resolution and annotated by top differentially expressed genes — an approach that is biased and resolution-dependent. CLANP provides an annotation-independent way to check whether a cluster label is biologically meaningful, by:

1. Running cNMF across a sweep of resolutions (K).
2. Scoring each gene expression program (GEP) against user-supplied cell-type labels via GSEA, producing a GEP × cell-type NES matrix.
3. Linking GEPs across K into a graph using KL divergence on NES vectors — programs that track the same biology persist across resolutions.
4. Labeling each GEP as a specific cell type, a family of related types, or ambiguous.
5. Scoring purity with an entropy-based score (ROGUE) on both the original labels and the GEP-derived groups.
6. Recommending an action per label: **keep** (already pure), **split** (label hides substructure), or **merge** (two labels are one program).

The output, per GEP at a given K: a cell-type or family label, a purity score, a recommended action, and the K range that resolves it.

## Why use it

- **Annotation-independent purity check.** ROGUE-based purity is computed from expression entropy, not marker genes.
- **Resolution-aware.** Persistence across K reveals which cell types are robustly resolvable and at what K.
- **Separates identity from activity.** GEPs distinguish what a cell *is* from what it is *doing*.
- **Actionable.** Each label gets a keep / split / merge recommendation grounded in program-level evidence.

## Installation

```bash
git clone https://github.com/dinh-lab/dinh-lab-clanp.git
cd dinh-lab-clanp
conda env create -f environment.yml
conda activate clanp
pip install -e .
```

Requires Python ≥ 3.9. Core dependencies: `scanpy`, `cnmf`, `gseapy`, `rogue` (or the Python port), `networkx`, `numpy`, `pandas`, `scipy`.

## Quick start

```python
import scanpy as sc
from clanp import run_pipeline

adata = sc.read_h5ad("myeloid_atlas.h5ad")

results = run_pipeline(
    adata,
    label_key="cell_type",          # column in adata.obs with user labels
    family_map="families.yaml",     # cell-type → family pools (optional)
    k_range=range(3, 36),           # resolutions to sweep
    out_dir="clanp_out/",
)

results.summary()        # per-label keep / split / merge calls
results.persistence()    # K-range each label is resolved over
results.gep_tree()       # hierarchical GEP graph across K
```

## Pipeline

```
scRNA-seq count matrix
        │
        ▼
1. cNMF across K (3–35)        → per-K GEP usage matrix
        │
        ▼
2. Score GEPs vs. labels (GSEA) → NES matrix (GEPs × cell types)
        │
        ▼
3. Build GEP graph across K     → edges weighted by KL divergence of NES
        │
        ▼
4. Label each GEP               → specific / family / ambiguous
        │
        ▼
5. Read persistence             → which GEPs are stable across K
        │
        ▼
6. ROGUE purity                 → per label and per GEP group
        │
        ▼
7. Compare to baseline          → label vs. best-matching GEP (Jaccard, ΔROGUE)
        │
        ▼
Output: per-label action — keep / split / merge, with K range
```

See `docs/pipeline.md` for the full method description.

## Inputs

- **Count matrix** — `AnnData` with raw counts in `adata.X` or `adata.raw.X`.
- **Cell-type labels** — a column in `adata.obs` (any string labels).
- **Family map** *(optional)* — YAML mapping cell types to higher-level families, used for the family-level fallback when a GEP is shared across related types. Example:

  ```yaml
  Monocyte_family: [Mono_CD14, Mono_CD16, Mono_CD14_IL1B, Mono_CD14_ID1, ...]
  Macrophage_family: [Mac_SPP1, Mac_IL1B, Mac_CXCL9, Mac_IL1Bint]
  DC_family: [cDC1_CLEC9A, cDC2_CD1C, cDC2_CD33, mregDC_LAMP3]
  ```

## Outputs

Written to `out_dir/`:

- `gep_labels.tsv` — per (K, GEP) row with assigned label, family, and ambiguity flag.
- `rogue_scores.tsv` — purity for each original label, each GEP group, and the label∩GEP intersection.
- `persistence.tsv` — first_K, last_K, and captured K values for each label.
- `recommendations.tsv` — per-label keep / split / merge call with supporting evidence.
- `figures/` — persistence plot, GEP tree, per-label diagnostic panels.

## Reproducing the paper

The head-and-neck cancer myeloid analysis (26,444 cells, 63 patients, 17 cell types) from the paper is reproducible from:

```bash
bash scripts/run_hnc_myeloid.sh
```

Inputs and expected outputs are described in `docs/reproduce.md`.

## Citation

If you use CLANP in your work, please cite:

> Sanjeev S., Dinh H.Q. *CLANP: Cluster Label Assessment via NMF Programs for scRNA-seq.* (in preparation).


## Contact

Issues and pull requests welcome. For research questions, contact Sumedha Sanjeev or Huy Q. Dinh at the McArdle Laboratory for Cancer Research, UW–Madison.
