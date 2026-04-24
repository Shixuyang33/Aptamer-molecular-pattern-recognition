# Reproducible analysis scripts for scApt-seq figure generation

This repository contains the analysis scripts associated with the manuscript **"Molecular Pattern Recognition by Sequencing Multiplex Aptamers for Cell Phenotyping at Single-cell Resolution"**. The scripts were organized to reproduce the figure-level analyses for the three major single-cell datasets used in the study:

* **Figure 4 / related supplementary figures**: HCT-8 cell-line validation dataset
* **Figure 5 / related supplementary figures**: ascites single-cell dataset
* **Figure 7 / related supplementary figures**: drug-resistance reversal dataset

These scripts are intended to provide a figure-aligned, manuscript-consistent analysis workflow that can be version-controlled and shared through GitHub.

## Repository contents

### 1. `Figure4_cellline_scAptseq_GitHub.R`

Analysis workflow for the **cell-line validation dataset** corresponding to **Figure 4** and related supplementary plots.

This script generates analyses corresponding to:

* **Figure 4a**: single-cell quality-control metrics across the four HCT-8 test groups
* **Figure 4b**: transcriptome UMAP, group identity, *PTK7* RNA, sgc8c aptamer signal, and scrambled aptamer control
* **Figure 4c**: relative proportions of nine target-specific aptamers plus the scrambled control across the four test groups
* **Figure 4d-e**: PTK7 target-dependent validation using the knockout/control dataset
* **Supplementary plots** related to HCT-8 aptamer-target feature visualization and cluster-level RNA-aptamer concordance

### 2. `Figure5_ascites_scAptseq_Fig5_Supplementary_GitHub.R`

Analysis workflow for the **ascites dataset** corresponding to **Figure 5** and related supplementary plots.

This script generates analyses corresponding to:

* **Figure 5a**: transcriptome-based clustering and sample identity in ascites samples
* **Figure 5b**: mean aptamer-expression heatmap between P1 and P2
* **Figure 5c-d**: mRNA/aptamer UMAP projections and correlation analyses for selected representative targets
* **Figure 5e-f**: patient-level and cluster-level violin plots for representative aptamers
* **Figure S7**: QC metrics, marker-gene dot plot, aptamer-based clustering, and patient separation
* **Figure S8**: additional cluster-level RNA-aptamer correlation analyses

### 3. `Figure7_drug_resistance_scAptseq_GitHub.R`

Analysis workflow for the **drug-resistance reversal model** corresponding to **Figure 7** and related supplementary plots.

This script generates analyses corresponding to:

* **Figure 7a**: transcriptome-based clustering across HCT-8, 5-FU, and 5-FU-Evo groups
* **Figure 7b**: aptamer-expression-based clustering across the three groups
* **Figure 7c**: heatmap of mean aptamer expression across the three groups
* **Figure 7d-f**: RNA, aptamer, and bulk-protein comparisons for representative targets
* **Figure S10-S14**: QC metrics, feature plots, RNA-aptamer correlation plots, transcript boxplots, and aptamer boxplots

## General analysis framework

Across the three scripts, the workflows follow the same general structure:

1. Load transcriptome and aptamer matrices
2. Perform quality control and filtering
3. Normalize and scale data
4. Integrate groups when required
5. Perform dimensional reduction and clustering
6. Generate figure-aligned visualizations and summary tables
7. Save intermediate Seurat objects and final figure files

## Expected inputs

The scripts assume access to processed single-cell transcriptome and aptamer count matrices generated from the scApt-seq workflow. Depending on the figure, the scripts may require:

* RNA count matrices
* aptamer count matrices
* optional metadata tables
* optional bulk proteomics table for Figure 7d-f
* optional PTK7 knockout/control dataset for Figure 4d-e

Input paths are controlled through **environment variables** rather than hard-coded absolute paths, so the scripts can be adapted across local workstations, servers, or cloud environments.


## Important notes

### 1. Figure-level consistency

These scripts were reorganized specifically to match the final manuscript structure. They are not meant to represent a universal end-to-end pipeline for all scApt-seq datasets; instead, they are tailored to reproduce the figure panels and related supplementary analyses described in the paper.

### 2. Aptamer naming and target mapping

The scripts include standardized aptamer naming and aptamer-target mappings to keep the analysis consistent with the manuscript text and figure legends. If your raw feature names differ from the current naming convention, you may need to update the feature-matching section of each script.

### 3. Sample naming

Please confirm that your local sample names match the mappings used in the scripts. In particular:

* ascites scripts assume the manuscript-level mapping for **P1** and **P2**
* drug-resistance scripts assume the three-group structure **HCT-8**, **5-FU**, and **5-FU-Evo**
* Figure 4 scripts assume the four HCT-8 test groups and, where applicable, the PTK7 control/knockout comparison

### 4. Proteomics comparison in Figure 7

The Figure 7 script treats bulk proteomics as an **orthogonal reference** rather than a single-cell readout. Protein values shown in Figure 7d-f are expected to come from **bulk TMT/LC-MS/MS** data from matched populations, not from the same single cells used in scApt-seq.

### 5. Correlation interpretation

Where RNA-aptamer correlations are computed, these values should be interpreted cautiously and in the biological context discussed in the manuscript. They are not intended to imply that aptamer-derived counts are universal quantitative surrogates of target abundance.

### 6. Environment and package versions

For reproducibility, it is recommended to record:

* R version
* package versions
* operating system
* session information

Each script is designed to save a `sessionInfo.txt` or related run summary where possible.

## Recommended package dependencies

Typical packages required include:

* `Seurat`
* `ggplot2`
* `dplyr`
* `patchwork`
* `pheatmap`
* `Matrix`
* `data.table`
* `stringr`
* `cowplot`
* `RColorBrewer`
* `harmony` (if integration is used)

Please install missing packages before execution.

## Reproducibility and GitHub use

To make the repository easier for others to use, it is recommended to include:

* this `README.md`
* the three figure-level R scripts
* a `data/` directory description or placeholder
* an `environment.yml` or package list
* a short note describing which raw or processed matrices are required but not publicly redistributed here

## Citation

If these scripts are used, please cite the associated manuscript.

## Disclaimer

These scripts were curated to align with the final submitted manuscript figures. Small adjustments may still be required if local matrix formats, barcode conventions, or metadata column names differ from those used during manuscript preparation.
