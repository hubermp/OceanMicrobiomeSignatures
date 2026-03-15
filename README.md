# OceanMicrobiomeSignatures

Analysis code associated with:

> **Huber et al.** – *Ocean Microbiome Signatures capturing the spatial organization of metabolic potential in the global sunlit ocean*  (under review)

---

## Repository Structure

| File | Description |
|------|-------------|
| `WF01_OPUs_final.ipynb` | Workflow 1: OPU definition, biogeography, and functional diversity (Python) |
| `WF02_OMS_final.ipynb` | Workflow 2: OMS identification, characterization, and mapping (Python) |
| `SDM_OMS.R` | Species Distribution Models for each OMS (R) |

---

## Input Data (Required)

All input tables are available as Supplementary Tables in the manuscript or on Zenodo (DOI: *to be assigned upon publication*). Place them in the same directory as the notebooks/scripts.

| Filename | Description | Supp. Table |
|----------|-------------|-------------|
| `TableS1_metadata.csv` | Metadata for 1,379 metagenomic samples | Table S1 |
| `TableS2_OPUs_abundance.csv` | OPU abundance matrix (419,979 × 1,379), normalized and rarefied | Table S2 |
| `TableS3_OPUs_annotation.csv` | OPU taxonomic (GTDB) and functional (KEGG KO) annotation | Table S3 |
| `TableS5_KeyKOdb.csv` | KeyKOdb – marker KEGG Orthology reference database | Table S5 |
| `TableS7_MAGs.txt` | 20,482 marine MAGs used for genome size analysis | Table S7 |
| `TableS8_OMS_W_matrix.csv` | NMF W matrix – OPU weights per OMS (419,979 × 10) | Table S8 |
| `H_10.tsv` | NMF H matrix – OMS contributions per sample (10 × 1,379) | Table S9 |
| `W_10.tsv` | NMF W matrix – raw TSV format | Table S8 |
| `Pres_Abs.txt` | Presence/absence matrix per OMS per sample | Table S11 |
| `Variables_ready/` | Environmental raster layers (.tif) from World Ocean Atlas 2018 | — |
| `M02.tif` | Mixed Layer Depth raster (WOA18) | — |

---

## Workflow 1 – `WF01_OPUs_final.ipynb`

Covers the full **Operational Protein Unit (OPU)** analysis pipeline.

**Sections:**
1. Load required packages
2. Load input tables (Tables S1–S5, S7)
3. Spatial distribution of metagenomes (Supplementary Fig. S1)
4. OPU functional characterization and global distribution
   - 4.1 Known vs. Unknown (dark matter) OPUs
   - 4.2 Global distribution: distance decay, ocean basin composition, Core vs. Endemic OPUs (Figs. 2, Table S4, S6)
   - 4.3 Functional diversity and richness across latitude (Supplementary Fig. S7)

**Note:** Mantel correlogram (Section 4.2.1) and UpSet plot (Section 4.2.3) were run in R; code is provided as markdown cells within the notebook.

---

## Workflow 2 – `WF02_OMS_final.ipynb`

Covers the identification and characterization of **Ocean Microbial Signatures (OMSs)** using Non-negative Matrix Factorization (NMF).

**Sections:**
1. Load H and W matrices from NMF (Tables S8–S9)
2. Normalize H matrix and merge with metadata
3. Geographic distribution of OMSs (major-contributing and monodominant; Table S11, Fig. 3)
4. Functional characterization (KeyKO, PCA, Z-score analysis; Tables S12–S13, Fig. 4)
5. Taxonomic characterization of OMSs (Fig. 5)
6. Environmental niche analysis (Table S10, Supplementary Fig. S4)

**Note:** NMF decomposition was performed externally using the custom OPUs pipeline ([github.com/celiosantosjr/OPUs_pipe](https://github.com/celiosantosjr/OPUs_pipe)). This notebook starts from the precomputed W and H matrices.

---

## SDM Script – `SDM_OMS.R`

Fits ensemble **Species Distribution Models (SDMs)** for each OMS to predict global Habitat Suitability Index (HSI) maps.

**Algorithms:** GLM, Random Forest (RF), ANN, Boosted Regression Trees (GBM/BRT), implemented via the `h2o` R package.  
**Evaluation:** 5-fold cross-validation; only models with AUC > 0.7 are retained in the ensemble.  
**Output:** Mean HSI rasters per OMS (`.tif`), binarized maps (TSS threshold), variable importance, and AUC scores.

Corresponds to the **Predictive Models** section of the Methods and **Fig. 3** (HSI maps).

---

## Software Requirements

### Python (Workflows 1 & 2)
```
Python >= 3.8
pandas · numpy · scipy · seaborn · matplotlib
scikit-learn · scikit-posthocs · scikit-bio · pycirclize
```
```bash
pip install pandas numpy scipy seaborn matplotlib scikit-learn scikit-posthocs scikit-bio pycirclize
```

### R (SDMs and auxiliary analyses)
```r
install.packages(c("raster", "terra", "usdm", "sf", "dplyr", "MASS",
                   "h2o", "dismo", "reshape", "ggplot2", "tidyr", "maps",
                   "tidyterra", "ggpubr", "gridExtra", "cowplot",
                   "RColorBrewer", "GGally", "viridis", "vegan",
                   "picante", "UpSetR"))
```

---

## Data Availability

Raw metagenomic data were processed using **MGnify v5.0** ([ebi.ac.uk/metagenomics](https://www.ebi.ac.uk/metagenomics/)).  
The full dataset is compiled in **AtlantECO-BASEv1** ([zenodo.org/records/7944433](https://zenodo.org/records/7944433)).  
Protein clustering used a custom pipeline: [github.com/celiosantosjr/OPUs_pipe](https://github.com/celiosantosjr/OPUs_pipe).  
Environmental layers: **World Ocean Atlas 2018** ([ncei.noaa.gov/products/world-ocean-atlas](https://www.ncei.noaa.gov/products/world-ocean-atlas)).

---

## Citation

Huber et al. (2026). *Ocean Microbiome Signatures capturing the spatial organization of metabolic potential in the global sunlit ocean*. *Nature Microbiology*. Under review.
