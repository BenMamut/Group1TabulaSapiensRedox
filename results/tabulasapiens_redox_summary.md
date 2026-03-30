# Tabula Sapiens Redox Transcriptional Programs — Complete Project Summary

## Research Question

**Are redox transcriptional programs primarily defined by cell type or donor-specific variation, and are distinct redox modules regulated independently or coordinately?**

Sub-questions:
1. Are redox transcriptional programs primarily intrinsic to cell type, or do they exhibit measurable donor-specific variation within the same cell type?
2. Do distinct redox modules (ROS production, detoxification, NADPH generation, signaling, metabolic coupling) vary independently, or are they tightly coordinated?
3. If a donor exhibits elevated transcription of a given redox module in one tissue or cell type, is this pattern coherent across other tissues?
4. How are redox modules coupled to other major biological programs (inflammation, proliferation, metabolism)?
5. Do redox-related transcriptional programs show more or less donor-level variability than other functional gene sets?

## Critical Limitation

These scores are **mRNA transcript-based proxies**, not direct biochemical measurements. Any redox score is separated from true redox potential by at least two regulatory layers:

> Transcript → Enzyme Activity → Metabolite Ratio → Redox Potential

We cannot conclude that cells are functionally more or less oxidized/stressed — only that **transcriptional programs associated with redox are more or less active**.

---

## Data Source

- **Database**: Tabula Sapiens via CELLxGENE Census (version 2025-11-08)
- **Access**: Direct S3 URI (`s3://cellxgene-census-public-us-west-2/cell-census/2025-11-08/soma/`)
- **Collection ID**: `e5f58829-1a66-40b5-a624-9046778e74f5`
- **Filters**: Primary-data cells only, Tabula Sapiens collection only
- **Donors**: 24 individuals
- **Tissues**: 30 tissue groups
- **Cells**: ~1.1 million total (uncapped run)

## Gene Modules (12 modules, 167 unique genes)

| Category | Module | # Genes | Example Genes |
|----------|--------|---------|---------------|
| **I. ROS Production** | NOX_SYSTEM | 17 | NOX1-5, DUOX1/2, CYBA, NCF1/2/4 |
| | MITO_ROS | 20 | NDUFA1/2/9, NDUFS1-3, SDHA-D, UQCRC1/2 |
| **II. ROS Detoxification** | SOD_SYSTEM | 4 | SOD1, SOD2, SOD3, CCS |
| | CATALASE_PEROXISOME | 18 | CAT, PRDX5, ACOX1-3, PEX family |
| | GLUTATHIONE_SYSTEM | 23 | GCLC, GCLM, GPX1-8, GST family, GLRX |
| | THIOREDOXIN_SYSTEM | 17 | TXN, TXNRD1-3, PRDX1-6, SRXN1 |
| **III. Reducing Power** | PPP_MODULE | 10 | G6PD, PGD, TKT, TALDO1 |
| | NADPH_DEHYDROGENASES | 12 | ME1-3, IDH1/2, MTHFD family, NNT |
| **IV. Redox Signaling** | NRF2_TARGETS | 19 | NFE2L2, KEAP1, NQO1, HMOX1 |
| | FOXO_OXSTRESS | 17 | FOXO1/3/4, SESN1-3, GADD45A/B |
| **V. Metabolic Coupling** | OXPHOS_MODULE | 16 | Complex I-V subunits, CYCS, TFAM |
| | GLYCOLYSIS_MODULE | 19 | HK1/2, PFKL/M/P, GAPDH, PKM, LDHA/B |

## Processing Pipeline

1. **Per-cell normalization**: CPM-like scaling (target sum 10,000) → log1p transform
2. **Per-cell module scoring**: Mean of log-normalized expression values for genes in each module
3. **Pseudobulk aggregation**: Group cells by donor_id × cell_type × tissue_general, compute mean scores per group (minimum 20 cells per group)

---

## Results Directory Guide

The `results/` folder contains outputs from multiple stages of analysis. Here is what each folder/file contains, whether it represents a **preliminary test** or a **final result**, and whether it should be used going forward.

### Preliminary / Superseded (can be ignored for final analysis)

These folders were from early development, iterative testing, or capped runs that have been superseded by the full uncapped analysis.

| Folder/File | What It Was | Status |
|-------------|-------------|--------|
| `preview_*.csv` | Initial Census queries testing which donors, cell types, tissues existed | **Preliminary** — exploratory metadata checks |
| `blood_BR1051_scores_2000.parquet` | Single-donor blood test (2000 cells) | **Preliminary** — pipeline proof-of-concept |
| `blood_3donors_*.parquet` | 3-donor blood test runs | **Preliminary** — early scaling tests |
| `analysis_BR1051_2000/` | Correlation/variance analysis on single-donor blood | **Preliminary** — single-donor, uninformative for donor effects |
| `analysis_blood_3donors_4000/` | Correlation/variance on 3-donor blood | **Preliminary** — too few donors for meaningful inference |
| `by_tissue_3donors_1500/` | Per-tissue parquets for blood/liver/lung with 3 donors | **Preliminary** — early multi-tissue test |
| `dryrun_combined_3tissues.parquet` | Combined pseudobulk from 3-tissue dry run (90 rows) | **Preliminary** — proof that combine step works |
| `dryrun_variance_3tissues/` | Variance partition on 3-tissue dry run | **Preliminary** — useful to confirm method, but insufficient scope |
| `tabulasapiens_batch_pseudobulk_3000/` | Full 30-tissue batch with 3000-cell cap per tissue | **Superseded** — replaced by uncapped run |
| `tabulasapiens_combined_pseudobulk_3000.parquet` | Combined capped table (649 rows, 23 donors, 120 cell types) | **Superseded** |
| `tabulasapiens_pseudobulk_combined_3000.parquet` | Duplicate/earlier version of combined capped table | **Superseded** |
| `tabulasapiens_variance_partition_3000/` | Variance partition on capped data | **Superseded** by full run |
| `variance_partition_tabulasapiens_3000/` | Empty — failed earlier attempt | **Delete candidate** |
| `tabulasapiens_combined_analysis_3000/` | Module correlations + donor variability on capped data | **Superseded** by full run |
| `tabulasapiens_blood_score_pseudobulk_5000.parquet` | Blood-only capped test | **Preliminary** |
| `tabulasapiens_blood_score_pseudobulk_full.parquet` | Blood-only uncapped QC test | **Preliminary** — used to verify uncapped extraction before full run |
| `redox_quadrant_plots/` | First attempt at ROS vs Antioxidant scatter (manually split into 2 vs 4 modules) | **Superseded** — biologically flawed axis construction (see PCA section below) |

### Active / Final Results

These are the outputs that should be used for analysis and manuscript.

| Folder/File | What It Contains | Status |
|-------------|-----------------|--------|
| **`tabulasapiens_batch_pseudobulk_full/`** | Per-tissue pseudobulk parquets for all 30 tissues (uncapped, ~1.1M cells) | **Primary data** |
| **`tabulasapiens_combined_pseudobulk_full.parquet`** | Master combined table: 1662 rows, 24 donors, 161 cell types, 30 tissues | **Primary dataset for all downstream analysis** |
| **`tabulasapiens_combined_analysis_full/`** | Module correlations + donor variability ranking (full uncapped data) | **Final** |
| **`tabulasapiens_variance_partition_full/`** | Variance partition: cell_type vs tissue vs donor (full data) | **Final** |
| **`donor_coherence_full/`** | Cross-tissue donor coherence analysis | **Final** |
| **`scoring_method_experiment_blood_12000/`** | Mean vs PCA1 scoring comparison (blood, 12k cells) | **Final** — sensitivity analysis |
| **`scoring_benchmark_blood_12000/`** | Mean vs AUCell-style vs background-adjusted scoring (blood, 12k cells) | **Final** — sensitivity analysis |
| **`redox_pca_plots/`** | PCA scatter plots using all 12 modules, raw and donor-corrected | **Final** — replaces redox_quadrant_plots |

---

## Key Findings

### 1. Variance Partition (addresses Sub-question 1)

**Cell type is the dominant driver of redox transcriptional variation.**

Full uncapped results (1662 rows, all 30 tissues):

| Factor | Mean Variance Fraction | Interpretation |
|--------|----------------------|----------------|
| Cell type | **64.3%** | Overwhelmingly dominant |
| Tissue | 5.3% | Moderate/contextual |
| Donor | 3.9% | Small but detectable |
| Residual | 26.5% | Noise + unmodeled factors |

Per-module breakdown:

| Module | Cell Type | Tissue | Donor | Residual |
|--------|-----------|--------|-------|----------|
| NOX_SYSTEM | 77.0% | 5.3% | 2.5% | 15.2% |
| GLUTATHIONE_SYSTEM | 75.7% | 2.6% | 2.7% | 19.0% |
| THIOREDOXIN_SYSTEM | 70.3% | 3.6% | 2.6% | 23.5% |
| NRF2_TARGETS | 66.4% | 3.2% | 5.2% | 25.1% |
| OXPHOS_MODULE | 66.4% | 4.7% | 2.4% | 26.5% |
| MITO_ROS | 64.6% | 5.7% | 2.7% | 27.0% |
| SOD_SYSTEM | 62.4% | 5.5% | 3.1% | 29.0% |
| CATALASE_PEROXISOME | 60.7% | 9.0% | 3.2% | 27.1% |
| PPP_MODULE | 61.1% | 8.2% | 3.8% | 26.8% |
| NADPH_DEHYDROGENASES | 64.3% | 4.3% | 3.8% | 27.5% |
| FOXO_OXSTRESS | **49.4%** | 6.1% | **9.2%** | 35.3% |
| GLYCOLYSIS_MODULE | 53.4% | 5.0% | 5.4% | 36.3% |

Notable: **FOXO_OXSTRESS** has the highest donor fraction (9.2%) and lowest cell-type fraction (49.4%) — this stress-response pathway shows the most individual-level variation. Glycolysis is the second most donor-variable module.

### 2. Module Coordination (addresses Sub-question 2)

**Most modules are positively correlated. They do NOT vary independently.**

Key Spearman correlations (full uncapped data):

| Module Pair | ρ | Interpretation |
|-------------|---|---------------|
| MITO_ROS ↔ OXPHOS | **0.987** | Nearly identical — measuring the same mitochondrial program |
| GLUTATHIONE ↔ THIOREDOXIN | **0.881** | Tightly coupled thiol-based defense systems |
| GLUTATHIONE ↔ NRF2_TARGETS | **0.848** | NRF2 transcriptionally drives glutathione genes |
| CATALASE ↔ MITO_ROS | 0.811 | Active mitochondria co-express peroxisomal defense |
| CATALASE ↔ NADPH_DEHYDROGENASES | 0.839 | Peroxisomal enzymes need NADPH supply |
| NOX_SYSTEM ↔ PPP_MODULE | 0.572 | NOX needs NADPH from PPP to produce ROS |
| NOX_SYSTEM ↔ SOD_SYSTEM | **-0.075** | Effectively uncorrelated — genuinely distinct axes |
| NOX_SYSTEM ↔ THIOREDOXIN | **-0.000** | Uncorrelated — innate immune ROS vs thiol defense |

**Biological interpretation**: There is a dominant "general metabolic activity" axis where cells that are metabolically active upregulate most redox modules simultaneously. NOX_SYSTEM is the most independent axis — it represents a dedicated innate-immune ROS production program that is largely uncoupled from the antioxidant/metabolic programs.

### 3. Cross-Tissue Donor Coherence (addresses Sub-question 3)

**Donor effects are modestly coherent across tissues** — a donor who is high for a module in one tissue tends to be somewhat high in others, but the effect is weak.

| Module | Mean Cross-Tissue Spearman | Interpretation |
|--------|---------------------------|----------------|
| GLYCOLYSIS | 0.259 | Most coherent — possibly systemic metabolic phenotype |
| GLUTATHIONE | 0.208 | |
| OXPHOS | 0.203 | |
| PPP | 0.199 | |
| NRF2_TARGETS | 0.195 | |
| CATALASE_PEROXISOME | 0.189 | |
| MITO_ROS | 0.185 | |
| NOX_SYSTEM | 0.184 | |
| FOXO_OXSTRESS | 0.161 | |
| NADPH_DEHYDROGENASES | 0.145 | |
| THIOREDOXIN | 0.125 | |
| SOD_SYSTEM | **0.076** | Least coherent — most tissue-specific donor effects |

These are modest but consistently positive correlations (237 tissue-pairs tested per module, requiring ≥4 shared donors per pair). The overall picture: donor-level redox variation is partially systemic, but mostly local to tissue context.

### 4. PCA of Redox Modules — What PC1 and PC2 Mean

Because the 12 modules are largely correlated (not independent), PCA provides the most informative 2D projection.

**PCA eigenvalue decomposition (all 12 modules, z-scored):**

| PC | Variance Explained | Biological Meaning |
|----|-------------------|-------------------|
| PC1 | **63.7%** | General redox metabolic intensity |
| PC2 | **11.6%** | NOX/innate-immune axis |
| PC3 | 8.5% | SOD-specific axis |
| PC4 | 4.7% | Minor variation |

**What PC1 represents ("Redox Metabolic Intensity"):**
All 12 modules load positively on PC1 with similar weights (+0.11 to +0.33). PC1 captures whether a cell has high or low overall transcription of redox-related genes. This is the dominant axis of variation: metabolically active cells (macrophages, monocytes, epithelial cells) score high on PC1; quiescent cells (T cells, erythrocytes) score low. **PC1 does NOT separate ROS production from antioxidant defense** — it is the shared metabolic program where active cells upregulate everything together.

**What PC2 represents ("NOX/Innate Immune Axis"):**
PC2 is dominated by NOX_SYSTEM (+0.77 loading) and PPP_MODULE (+0.41), with antioxidant modules loading negatively (SOD -0.28, THIOREDOXIN -0.25). PC2 separates cells that rely on dedicated NADPH oxidase-based ROS production (professional phagocytes like neutrophils, monocytes) from cells that are more oriented toward antioxidant defense and mitochondrial metabolism. **PC2 is the closest axis to a "ROS production vs defense" separation**, but it specifically captures NOX-dependent innate immune ROS, not generic ROS.

**How to interpret the PCA scatter plot quadrants:**

| Position | PC1 | PC2 | Meaning | Cell Types |
|----------|-----|-----|---------|------------|
| Upper-right | High | High | High overall redox + strong NOX | Monocytes, macrophages |
| Lower-right | High | Low | High overall redox + strong antioxidant/mito | Fibroblasts, endothelial cells |
| Upper-left | Low | High | Low overall redox but NOX-active | Some immune subtypes |
| Lower-left | Low | Low | Low redox activity overall | T cells, erythrocytes |

**Why the original 2-axis (ROS vs Antioxidant) plot was flawed:**
The first attempt manually assigned 2 modules to "ROS production" and 4 to "Antioxidant defense" (ignoring the other 6). This was biologically misleading because:
- The 6 omitted modules all correlated positively with both axes (r = 0.47–0.81)
- The dominant signal (PC1 @ 63.7%) loads on ALL modules equally — cells don't separately regulate "ROS" and "antioxidant" at the transcriptional level; they scale the entire redox program up or down together
- This caused diagonal streaks (both axes measured the same underlying metabolic activity), and most cell types appeared in "low antioxidant" quadrants because the 2-vs-4 module split was asymmetric

The PCA approach correctly lets the data define its own axes. The result: lineages separate clearly (myeloid vs lymphoid vs stromal vs endothelial), and the NOX-driven immune axis emerges naturally as PC2.

### 5. Donor Effects on Cluster Tightness

**Donor variation adds real but modest noise** to cell-type clusters in PCA space.

Donor correction (subtracting per-donor mean shifts) was applied and compared to raw PCA:

| Metric | Value |
|--------|-------|
| Mean spread reduction (top-30 cell types) | **0%** (negligible overall) |
| Mean spread reduction (all cell types) | 1.7% |
| Best individual improvements | Endothelial lymphatic (-14%), Basal cell (-11%), Monocyte (-10%) |

**Interpretation**: Donor effects exist but they are dwarfed by the biological differences between cell types and by tissue-context effects (same cell type in different tissues). The "streaks" in the raw PCA scatterplot are primarily caused by tissue-context variation (e.g., a macrophage in lung vs liver vs fat has different redox programs), not by donor noise. This is itself a biologically informative finding.

### 6. Scoring Method Sensitivity

Three alternative scoring methods were benchmarked on blood tissue (12,000 cells):

**Mean vs PCA1 scoring (per-module):**
- Most modules show near-perfect agreement (Spearman > 0.9)
- Exception: SOD_SYSTEM (Spearman = 0.13) — only 4 genes, PCA1 captures a different aspect than the mean
- Conclusion: Mean scoring is robust for modules with ≥10 genes

**Mean vs AUCell-style vs Background-adjusted:**

| Comparison | Mean Spearman |
|------------|--------------|
| AUCell vs Background-adjusted | 0.831 |
| Mean vs Background-adjusted | 0.626 |
| Mean vs AUCell | 0.533 |

**Effect on variance decomposition (blood only):**

| Method | Cell-type fraction | Donor fraction |
|--------|-------------------|----------------|
| Mean | 83.5% | 6.5% |
| Background-adjusted | 84.8% | 5.2% |
| AUCell-style | 86.5% | 4.2% |

AUCell-style scoring pushes more variance to cell type and reduces apparent donor signal. The baseline mean score preserves the most donor-level information, which is desirable for our research question about donor variation.

**Recommendation**: Keep mean score as primary method. Report AUCell-style and background-adjusted as sensitivity analyses to demonstrate robustness of main conclusions.

### 7. Donor Variability Ranking

Across the full dataset, modules ranked by inter-donor variance:

| Module | Donor Variance | Interpretation |
|--------|---------------|----------------|
| SOD_SYSTEM | 0.125 | Most variable across donors |
| FOXO_OXSTRESS | 0.111 | Stress response — donor-specific |
| OXPHOS_MODULE | 0.093 | |
| MITO_ROS | 0.078 | |
| THIOREDOXIN_SYSTEM | 0.072 | |
| GLYCOLYSIS_MODULE | 0.064 | |
| NRF2_TARGETS | 0.060 | |
| PPP_MODULE | 0.043 | |
| GLUTATHIONE_SYSTEM | 0.036 | |
| NOX_SYSTEM | 0.021 | |
| CATALASE_PEROXISOME | 0.020 | |
| NADPH_DEHYDROGENASES | 0.019 | Least variable across donors |

SOD_SYSTEM (4 genes: SOD1, SOD2, SOD3, CCS) is the most donor-variable module but also the least cross-tissue coherent — it behaves differently in different tissues and donors.

---

## Pipeline Scripts

| Script | Purpose |
|--------|---------|
| `redox_census_pipeline.py` | Main pipeline: Census querying, scoring, pseudobulk, batch extraction, variance partition, correlation analysis |
| `redox_advanced_analyses.py` | Cross-tissue donor coherence + multi-method scoring benchmark |
| `redox_scoring_method_experiment.py` | Mean vs PCA1 scoring comparison |
| `redox_quadrant_plot.py` | Original 2-axis scatter (superseded by PCA approach) |
| `redox_pca_scatter.py` | PCA-based scatter plots (raw, donor-corrected, lineage-colored, biplot) |

---

## Summary of Answers to Research Sub-Questions (So Far)

| Sub-Question | Answer from Current Data |
|-------------|------------------------|
| **1. Cell type vs donor?** | Cell type dominates (64.3%), donor is small (3.9%) but detectable. FOXO_OXSTRESS is the most donor-variable module (9.2%). |
| **2. Independent or coordinated?** | Mostly coordinated — 63.7% of variance is a shared "redox intensity" axis. NOX_SYSTEM is the most independent module. MITO_ROS and OXPHOS are nearly redundant (ρ=0.987). |
| **3. Cross-tissue donor coherence?** | Modest positive coherence (mean Spearman 0.08–0.26). Glycolysis most coherent, SOD least. Donor redox programs are partially systemic but mostly tissue-contextual. |
| **4. Coupling to other programs?** | Not yet tested (requires integrating non-redox gene sets). |
| **5. Redox vs other gene sets?** | Not yet tested (requires control/random pathway baselines). |

## Remaining Work

1. **Control pathway comparison**: Compare variance decomposition of redox modules against random gene sets or housekeeping pathways to determine if redox programs are unusually conserved or variable
2. **Residual coupling analysis**: Remove cell-type and tissue effects from module scores and compute correlations on residuals to isolate donor-level coordination patterns
3. **Coupling to other biological programs**: Integrate inflammatory, proliferative, and metabolic gene sets to test how redox modules co-vary with non-redox axes
4. **Stability / downsampling checks**: Resample to matched donor counts across tissues and verify that findings are robust
