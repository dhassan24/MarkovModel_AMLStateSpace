
**Author:** Danya Hassan  
**Institution:** Keck Graduate Institute (KGI), MS in Human Genetics and 
Genomics Data Analytics  
**Capstone Partner:** City of Hope — Dr. Russell Rockne Lab  
**Acknowledgements:** Thank you to the City of Hope Department of Mathematical 
Oncology and Beckman Research Institute for their endless support and guidance.
D.Frankhouser, D.O’Meally, J.Ambriz, Z.Chen, YH.Fu,
J.Irizarry, B.Zhang, YH.Kuo, G.Marcucci, B.Fortini, C.Espinoza, R.Rockne 

---

## Overview

This repository reproduces and extends the transcriptomic state-space model of
acute myeloid leukemia (AML) described in Rockne et al. (2021), applied to the
2018 AML mouse cohort. Peripheral blood mononuclear cell (PBMC) mRNA time-series
data from a CBFβ-SMMHC (CM) fusion-gene mouse model are used to:

1. Construct a low-dimensional transcriptomic state-space via SVD
2. Build a quasi-potential double-well energy landscape along the dominant
   disease-progression axis (PC2)
3. Simulate stochastic transcriptomic trajectories using a
   Metropolis–Hastings Monte Carlo algorithm
4. Predict individual mouse time-to-leukemia (first-passage time to State 3)
5. Validate predictions at both the per-mouse and population level
6. Identify gene expression programs associated with state transitions using
   tradeSeq trajectory modeling and GO enrichment analysis

---

## Repository Structure

```
/notebooks          Quarto (.qmd) analysis documents — rendered to HTML
/R                  Shared helper functions (if applicable)
/data               Input data files (see Data Files section below)
/outputs            Generated CSVs, figures, and animation outputs
```

---

## Notebooks

The four notebooks should be run in order. Each document is self-contained
(all parameters are re-declared at the top) but later notebooks depend on
output files written by earlier ones.

### 1. `aml2018_state_space_construction.qmd`

Constructs the AML transcriptomic state-space from the 2018 mRNA cohort.

**What it does:**
- Loads the 2018 `SummarizedExperiment` object and filters low-expression genes
- Log₂-transforms and mean-centers TPM values across Ctrl and CM samples
- Performs SVD on the centered expression matrix to extract principal components
- Flips and rotates the PC1–PC2 plane so that Ctrl trajectories are flat and
  the CM (leukemic) direction corresponds to negative PC2 values
- Projects all samples into the rotated space and removes outlier samples
  (mouse 3368, mouse 3336 week 8, mouse 3357 post-week-3)
- Extracts PC2 eigengenes (gene loadings) and saves them for downstream use
- Clusters CM samples into three k-means states along PC2:
  K1 (perturbed hematopoiesis), K2 (transition region), K3 (deep leukemia)
- Estimates critical points c₁, c₂, c₃ and constructs the quartic
  quasi-potential double-well curve
- Verifies the curvature at each critical point and plots the Boltzmann ratio
  across candidate c₂ values

**Key outputs written:**
| File | Description |
|---|---|
| `eigengenes_PC2_aml2018.csv` | PC2 gene loadings (unscaled and variance-scaled) for all expressed genes |
| `AML_2018_StateSpace_Mice.csv` | Per-sample state-space coordinates merged with full sample metadata |
| `final_aml_mrna_state_space.csv` | Filtered state-space coordinates (outlier samples removed) |
| `cm2018_sp_kmeansclusters.csv` | CM mouse samples with k-means cluster labels (K1/K2/K3) |
| `aml_qp_dw_k2max.csv` | Double-well quasi-potential curve using c₂ = K2 upper bound |
| `amlss_quasipotential_dwcurve.csv` | Double-well quasi-potential curve using c₂ = K2 centroid |
| `quasipotential_varying_c2.csv` | Family of quasi-potentials across the full c₂ candidate range |

---

### 2. `aml2018_mc_parameterization.qmd`

Calibrates and tunes the Monte Carlo simulation parameters before running
full predictions.

**What it does:**
- Rescales the quasi-potential so that the S1→S2 energy barrier equals
  0.99 A.U., matching the reference paper's physical regime (D/barrier ≈ 0.45)
- Defines the Metropolis–Hastings proposal function and first-passage time (FPT)
  runner with reflecting boundary conditions
- Calibrates `step_to_month`: anchors the mean simulated FPT (from random
  State 1 starts) to the mean observed leukemia time of the four calibration
  mice (3334, 3336, 3341, 3357)
- Runs an 18-combination step-size sweep (`step_sd` × `max_mult`) and selects
  the proposal parameters that minimise mean MAE at ≥ 80% State-3 reachability
- Runs a D sensitivity sweep (D = 0.25–0.75) with independent recalibration at
  each D value; validation mice (3370, 3338) are held out throughout
- Confirms that D = 0.45 (paper value) is near-optimal on this cohort and that
  MAE is stable across D = 0.40–0.55

**Calibration / validation split:**
| Role | Mouse IDs |
|---|---|
| Calibration | 3334, 3336, 3341, 3357 |
| Validation | 3370, 3338 |

**Selected parameters:** `step_sd = 15`, `max_mult = 1.5`, `D = 0.45`

---

### 3. `aml2018_mc_predictions.qmd`

Runs the calibrated MC simulation to produce leukemia onset predictions and
constructs the Markov transition model.

**What it does:**
- Re-declares all landscape and simulation objects for reproducibility
- Runs 700 FPT simulations per mouse starting from each mouse's first observed
  PC2 coordinate; reports predicted mean, median, mode, and SD
- Plots per-mouse FPT distributions (histogram + observed time overlay)
- Builds a filtered eigengene loading bar chart after removing confounding gene
  classes (Y-chromosome, ribosomal, mitochondrial, immunoglobulin, etc.)
- Constructs τ-lagged Markov transition matrices at τ = 40, 50, 60 to confirm
  lag-time stability
- Computes analytic mean first-passage times (MFPT) via the fundamental matrix
  N = (I − Q)⁻¹ for both population-level (random State 1 start) and
  per-mouse matrices
- Computes the stationary distribution of the three-state chain
- Validates population-level predictions: 1,000 FPT simulations from the c₁
  centroid, bootstrap test confirming the observed cohort mean is statistically
  consistent with the model (all 6 observed times fall within predicted IQR)
- Generates an animated GIF of a single MC trajectory on the rescaled landscape

**Key results:**
| Metric | Value |
|---|---|
| Population MFPT (S1 → S3) | ~6–7 months |
| Population MFPT (S2 → S3) | ~3 months |
| Observed cohort mean | 5.75 months |
| Bootstrap p-value | 0.73 (not significant — consistent with model) |
| All 6 observed times within predicted IQR | ✓ |

**Key outputs written:**
| File | Description |
|---|---|
| `AML_MC_C3pred_rescaledU_D045.gif` | Animated MC trajectory on the rescaled quasi-potential |

---

### 4. `aml2018_tradeseq_analysis.qmd`

Fits gene expression trajectories along the AML pseudotime axis and links
divergence clusters to biological pathways.

**What it does:**
- Filters to the top 250 PC2 eigengenes (highest absolute variance-scaled
  loading) from `eigengenes_PC2_aml2018.csv`
- Maps Ensembl IDs to gene symbols via biomaRt
- Fits negative binomial GAMs along pseudotime for two lineages (Ctrl, CM)
  using tradeSeq (k = 8 knots selected via `evaluateK`)
- Clusters genes by their CM − Ctrl divergence trajectory using hierarchical
  clustering (Ward's method, k = 6), capturing distinct temporal programs
- Plots per-cluster trajectories: spaghetti, mean ± SE, and faceted overview
- Constructs expressed and variance-filtered gene universes from the full
  mrna_2018 object for statistically valid ORA backgrounds
- Runs GO Biological Process ORA per cluster (clusterProfiler, BH correction)
- Highlights curated terms of AML relevance: mast cell degranulation (C1),
  angiogenesis (C2), epigenetic/chromatin regulation (C4), ribosomal
  transcription (C6)
- Runs GSEA on a genome-wide LFC ranking (CM vs Ctrl) and maps core-enrichment
  genes back to divergence clusters
- Generates STRING protein interaction networks for selected clusters

---

## Data Files

### Input (not included in repository — required to run)

| File | Description | Used in |
|---|---|---|
| `AML_mRNA_2018_1_real.rds` | `SummarizedExperiment` object with TPM, counts, and scaled counts assays; 33,789 genes × 139 samples | Notebooks 1, 2, 3, 4 |
| `metadata_mmu.txt` | Sample-level metadata for all mouse assays (semicolon-delimited) | Notebook 1 |

### Generated (written by notebooks, used downstream)

| File | Written by | Used in | Description |
|---|---|---|---|
| `eigengenes_PC2_aml2018.csv` | Notebook 1 | Notebooks 3, 4 | PC2 gene loadings — unscaled (`eigengene_PC2`) and variance-scaled (`eigengene_PC2_scaled`) for all expressed genes |
| `final_aml_mrna_state_space.csv` | Notebook 1 | Notebooks 2, 3, 4 | Per-sample PC2 coordinates after outlier filtering; includes `Time`, `mouse_id`, `group_name`, `AML_space.calc` |
| `AML_2018_StateSpace_Mice.csv` | Notebook 1 | Notebooks 2, 4 | State-space coordinates joined to full sample metadata from `mrna18_info` |
| `cm2018_sp_kmeansclusters.csv` | Notebook 1 | — | CM mouse samples annotated with K1/K2/K3 cluster labels |
| `aml_qp_dw_k2max.csv` | Notebook 1 | Notebooks 2, 3 | Quasi-potential curve (500 grid points); columns `AML_PC2_Location`, `U_vals`, `state`; c₂ = K2 upper bound |
| `amlss_quasipotential_dwcurve.csv` | Notebook 1 | — | Alternative quasi-potential with c₂ = K2 centroid |
| `quasipotential_varying_c2.csv` | Notebook 1 | — | Family of 10 quasi-potential curves across candidate c₂ values |
| `aml_qp_dw_k2max_rescaled.csv`* | Notebook 2 | — | Quasi-potential rescaled to 0.99 A.U. barrier (intermediate diagnostic) |
| `AML_MC_C3pred_rescaledU_D045.gif` | Notebook 3 | — | Animated MC trajectory visualization |

*Generated internally during parameterization; the rescaling is re-applied at
the top of Notebook 3 so no separate file is strictly required.

---

## Methods Summary

### State-Space Construction
SVD is applied to the log₂-transformed, mean-centered TPM matrix of Ctrl and CM
samples. The right singular vectors (columns of V) are gene loadings
(eigengenes); the left singular vectors (columns of U) are sample coordinates.
The PC1–PC2 plane is rotated to flatten Ctrl trajectories, and PC2 is signed so
that leukemia corresponds to negative values.

### Quasi-Potential
The quasi-potential is a quartic polynomial whose derivative factors as
α(x − c₁)(x − c₂)(x − c₃), giving two stable minima (c₁, c₃) and one
unstable maximum (c₂). c₁ and c₃ are estimated as k-means cluster centroids;
c₂ is set to the upper boundary of the transition cluster (K2). The potential
is rescaled so the c₁ → c₂ barrier equals 0.99 A.U. per the reference paper.

### Monte Carlo Simulation
A Metropolis–Hastings random walk on the discretized landscape simulates
stochastic transcriptomic dynamics. Proposals are Gaussian in PC2 coordinate
space (σ = 15, clipped to ±22.5 units) with reflecting boundaries.
Acceptance probability is min(1, exp(−ΔU/D)) with D = 0.45. First-passage
time to State 3 is recorded as the predicted leukemia onset.

### Markov Modeling
τ-lagged transition matrices are estimated from long MC trajectories and used
to compute analytic mean first-passage times via the fundamental matrix
N = (I − Q)⁻¹, where Q is the 2 × 2 submatrix of transient states.

---

## Dependencies

```r
# Bioconductor
BiocManager::install(c(
  "SummarizedExperiment", "tradeSeq", "SingleCellExperiment",
  "slingshot", "clusterProfiler", "org.Mm.eg.db",
  "enrichplot", "biomaRt"
))

# CRAN
install.packages(c(
  "ggplot2", "dplyr", "tidyr", "tidyverse", "tibble",
  "patchwork", "readr", "stringr", "ggrepel",
  "gganimate", "gifski", "STRINGdb", "tictoc"
))
```

R version 4.3+ recommended.

---

## Reference

Rockne, R.C. et al. (2021). *State-space and quasi-potential model of
transcriptomic dynamics underlying AML disease progression.*
Used as the methodological basis for this work.