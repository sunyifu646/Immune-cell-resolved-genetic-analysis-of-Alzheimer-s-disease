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

---

## Citation

If you use this code or adapt this workflow, please cite:

```
Sun Y, Zhong Z, Zhou W, Zhao L, Xie C.

An integrative single-cell genetic analysis prioritizes immune cell-specific genes and drug repurposing hypotheses in Alzheimer's disease.
```
