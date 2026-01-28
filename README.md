# CIRI (CRISPR Interference and Activation) Analysis Pipeline

This repository contains an R-based pipeline for processing and analyzing single-cell RNA-seq data from CRISPRa/i screens. It handles perturbation deconvolution, quality control, dimensionality reduction (Monocle3), and trajectory analysis to assess the impact of guides on cell differentiation.

## Prerequisites
The pipeline runs in a Dockerized environment to ensure reproducibility.

Base Image: rocker/r-ver:4.4.2 (Ubuntu 24.04 LTS).

Key Libraries: Monocle3, Seurat, Tidyverse, Bioconductor.

### Build Docker Image

```docker build -t ciri_pipeline . ```

## Pipeline Overview
The analysis follows this sequence:

- Deconvolution: Assign CRISPR perturbations to cells based on UMI distributions.

- Filtering: Quality control and protein-coding gene selection.

- Loading: Normalization and initial dimensionality reduction.

- Validation: Check knockdown/activation efficiency.

- Cluster Analysis: Perturbation enrichment in target clusters.

- Trajectory Analysis: Pseudotime dynamics and statistical testing.

- Signatures: Biological scoring.

## Usage
### 1. Perturbation Deconvolution

Assigns guide identities based on "fixed" (Cas9 modality) and "variable" (target) guide capture. This is the first and mandatory experiment of the analysis.

The front end function were created with Baryon (see Fairflow-BioinformaticsFramework/baryon-lang). 

Command:

```./perturbation_assignment/CIRI_assign.sh  --target_directory "/home/AB11_screening_CRISPR/analysis_11/."  --script_directory "/home/CIRI_bioinfo/perturbation_assignment/."   --guide_strategy 1   --matrix_filename "filtered_feature_bc_matrix.h5"   --threshold_a -1   --threshold_i -1```

Option, Required, Description, Example
--target_directory, Yes, The working directory containing your input files (guides.csv and the .h5 matrix).Note: Output files will be saved here.,/data/experiment_1


--script_directory,Yes,The directory where the CIRI_unified.R script is located.,/scripts/bioinfo
--guide_strategy,Yes,The analysis logic mode.1 = Single Guide (Ratio based)2 = Dual Guide (Sum/Rank based),2
--matrix_filename,Yes,The specific filename of the 10x Genomics H5 matrix inside the target directory.,filtered_feature_bc_matrix.h5
--threshold_a,No,Manual UMI threshold for CRISPRa. Set to -1 to auto-calculate.,23.5 (or -1)
--threshold_i,No,Manual UMI threshold for CRISPRi. Set to -1 to auto-calculate.,-1

Methodology:

Fixed Guides: Fixed guides are assessed to establish the perturbation type (CRISPRa vs CRISPRi). UMI counts are summarized into histograms (2000 bins for CRISPRa, 6000 for CRISPRi) and smoothed to identify peaks. The threshold is defined as the valley between the noise peak and the signal peak.

Variable Guides (Single): Variable guides are assigned if the top guide has ≥10 UMIs and the ratio between the top and second guide is ≥5.

Variable Guides (Dual): The sum of the top two guides must be ≥4 UMIs, and the ratio between this sum and the third guide must be ≥10. Cells are filtered if the two variable guides do not target the same gene.


### 2. Annotation & Filtering

Performs quality control on the gene expression matrix.

Filters Applied:

- Percentage of ribosomal and mitochondrial genes is calculated.

- Cells with <250 total called genes are filtered out.

- Mitochondrial and ribosomal genes are removed, keeping only protein-coding genes.

- Genes with <3 total UMIs are removed.

Command:

``` Rscript anno_filter.R /analysis/data/ filtered_feature_bc_matrix.h5 ```

### 3. Data Loading & Preprocessing

Initializes the Monocle3 object, performs normalization, and generates the initial UMAP.

Methodology:

Monocle3 clustering is performed. 

Command:

Args: <directory> <annotated_matrix_csv> <clustering_resolution>

``` Rscript CIRI_load.R /analysis/data/ annotated_matrix.csv 1e-5 ```

### 4. Target Validation

Validates perturbation efficiency by comparing target gene expression in perturbed cells vs. non-targeting controls.

Command:

``` Rscript guide_genes_expr.R ```

### 5. Cluster Enrichment Analysis

Analyzes the distribution of perturbations across clusters (specifically the target muscle cluster).

Methodology:

- The cluster of interest is selected.

- The percentage of cells in the target cluster is calculated for each guide combination.

- Combinations are filtered for a minimum cell count (e.g., 10 or 20).

- Combinations present in both samples with sufficient cells are visualized on heatmaps and scatterplots.

Command:

Args: <directory> <cluster_ids> <control_name> <min_cells>

``` Rscript CIRI_sec.R /analysis/data/ 4 "NTCa-NA" 40 ```

### 6. Subclustering & Trajectory

Subsets the target lineage (e.g., muscle cluster), re-clusters at high resolution, and learns pseudotime trajectories.

Methodology:

- The cluster of interest is subclustered.

- Pseudotime is calculated using Monocle3, setting the root at the node with the highest SOX2 expression.

Command:

Args: <directory> <cluster_ids> <root_gene> <group_name> <resolution>

``` Rscript CIRI_sec_subclusters.R /analysis/data/ 4 SOX2 muscle 1e-4 ```

### 7. Pseudotime Statistics (KS Test)

Performs Kolmogorov-Smirnov tests to detect significant shifts in differentiation speed (pseudotime distribution) compared to controls.

Command:

Args: <dir> <cds_rdata> <pseudotime_csv> <analysis_name> <run_all> <control_grp> <min_cells>

``` Rscript CIRI_pseudotime.R /analysis/data/ cds_sample_1_ordered.RData pseudotime_muscle_1.csv muscle_1 TRUE "NTCa_1A;NTCa_1B-NA" 8 ```

### 8. Signature Scoring

Scores cells based on predefined gene sets (e.g., Cell Cycle, Sarcomere Core) defined in signatures.R.

``` Rscript signatures.R ```


## Input/Output Structure
Input:

- filtered_feature_bc_matrix.h5 (10x Genomics output)

- guides.csv (Guide library definition)

Output:

- processed_cds.RData: Monocle3 object.

- annotation_data.csv: Cell-guide assignments.

- pseudotime_*.csv: Pseudotime values per cell.

- *.pdf: UMAPs, Heatmaps, Violin plots, and Volcano plots.
