# Immune-cell-resolved-genetic-analysis-of-Alzheimer-s-disease

## Overview

This repository contains the source code for the study:

**"An integrative single-cell genetic analysis prioritizes immune cell-specific genes and drug repurposing hypotheses in Alzheimer's disease"**

This project provides an integrative computational framework to investigate immune cell-specific genetic regulation associated with Alzheimer's disease (AD) susceptibility.

The workflow integrates multiple layers of genetic and transcriptomic evidence, including:

- single-cell expression quantitative trait loci (sc-eQTL) analysis
- two-sample Mendelian randomization (MR)
- Bayesian colocalization analysis
- independent immune cell eQTL replication
- peripheral blood single-cell RNA sequencing (scRNA-seq)
- drug signature enrichment analysis
- molecular docking analysis
- phenome-wide association study (PheWAS)

The framework aims to identify immune cell-resolved genetic regulators of Alzheimer's disease and generate exploratory drug repurposing hypotheses.

---

## Study Workflow

The analytical pipeline includes the following major steps:

```
Immune cell-specific sc-eQTL datasets
                |
                v
Two-sample Mendelian Randomization
                |
                v
Bayesian Colocalization Analysis
                |
                v
Independent DICE sc-eQTL Replication
                |
                v
Peripheral Blood scRNA-seq Characterization
                |
                v
Drug Signature Enrichment Analysis
                |
                v
Molecular Docking Analysis
                |
                v
Phenome-wide Association Study (PheWAS)
```

---

## Repository Structure

```
.
├── README.md
├── scripts/
│   ├── 01_MR_analysis/
│   ├── 02_colocalization/
│   ├── 03_DICE_replication/
│   ├── 04_scRNAseq_analysis/
│   ├── 05_drug_repurposing/
│   ├── 06_molecular_docking/
│   └── 07_PheWAS/
│
├── data/
│   └── data_access_information.md
│
├── results/
│   └── figures_and_tables/
│
└── environment/
    ├── R_packages.txt
    └── conda_environment.yml
```

---

## Data Availability

All analyses in this study were performed using publicly available datasets.

The source code supporting the findings of this study is available at:

https://github.com/sunyifu646/Immune-cell-resolved-genetic-analysis-of-Alzheimer-s-disease

The datasets used in this study can be obtained from their original public repositories:

### Immune cell sc-eQTL datasets

**OneK1K cohort**

- Immune cell-specific cis-eQTL summary statistics
- Used for discovery-stage MR analyses

**DICE project**

- Independent immune cell eQTL resource
- Used for replication analyses

### Alzheimer's disease GWAS

Publicly available Alzheimer's disease genome-wide association study summary statistics were used as outcome data.

### Peripheral blood scRNA-seq dataset

**GEO accession: GSE181279**

- Peripheral blood single-cell RNA sequencing dataset
- Used for exploratory characterization of prioritized genes

Detailed accession information and preprocessing procedures are described in the manuscript.

---

## Software Requirements

### R Environment

The main statistical analyses were performed using:

```
R version 4.4.2
```

Major R packages include:

- TwoSampleMR
- coloc
- data.table
- tidyverse
- ggplot2
- Seurat

A complete list of required R packages is provided in:

```
environment/R_packages.txt
```

---

## Reproducibility

To reproduce the analyses:

### Step 1. Download required datasets

Download all publicly available datasets from their original repositories.

### Step 2. Configure computational environment

Install required R packages:

```
environment/R_packages.txt
```

or configure the conda environment:

```
environment/conda_environment.yml
```

### Step 3. Run analysis pipelines

The scripts are organized according to the analytical workflow:

```
01_MR_analysis
        |
        v
02_colocalization
        |
        v
03_DICE_replication
        |
        v
04_scRNAseq_analysis
        |
        v
05_drug_repurposing
        |
        v
06_molecular_docking
        |
        v
07_PheWAS
```

The scripts generate intermediate files, statistical results, and figures corresponding to the analyses reported in the manuscript.

---

## Main Outputs

This repository enables generation of:

- immune cell-specific AD-associated eGenes
- Mendelian randomization results
- Bayesian colocalization results
- independent replication results
- single-cell expression characterization of prioritized genes
- drug–gene enrichment results
- molecular docking predictions
- PheWAS profiles for genetic safety assessment

---

## Key Findings Supported by This Code

The analytical framework identifies immune cell-specific genetic signals associated with Alzheimer's disease and prioritizes candidate genes including:

- KANSL1-AS1
- CTSH
- FCER1G
- EPHA1-AS1
- CRHR1
- JAZF1

The drug repurposing and molecular docking analyses provide computational hypotheses for future experimental validation.

---

## Citation

If you use this code or adapt this workflow, please cite:

```
Sun Y, Zhong Z, Zhou W, Zhao L, Xie C.

An integrative single-cell genetic analysis prioritizes immune cell-specific genes and drug repurposing hypotheses in Alzheimer's disease.
```

---

## License

This repository is provided for academic and non-commercial research purposes.

Please refer to the LICENSE file for details.

---
