# Molecular caracterisation of the sarcomatoid (SARC) bladder cancer variant 
These scripts were used to perform transcriptomic and genomic profiling of the sarcomatoid bladder cancer variant in paired samples. This project is conducted by the IMPACT team at Institut Curie (Paris, France). 

## Description
These scripts analyze 3'RNA-seq and WES data from the sarcomatoid bladder cancer variant project (IMPACT team, Institut Curie). The transcriptomic pipeline covers preprocessing, visualization, molecular class prediction, differential expression, GSEA, and microenvironment deconvolution.

**Note (2026-08-17):** The WES pipeline is not yet available on this repository and is still under development.

## Installation
```bash
git clone https://github.com/mgautheret/SARC_molecular_caracterisation.git
```

## Transcriptomic pipeline overview

Overall, the pipeline takes as input raw RNA-seq counts, TPM counts, and sample metadata, and produces classified/cleaned tables, differential expression results, and various QC/analysis plots.

| Script | Description |
|---|---|
| **MolecularClassif_NMIBC+MIBC.Rmd** | Classifies samples into MIBC/NMIBC subtypes using the [consensusMIBC](https://github.com/cit-bioinfo/consensusMIBC) and [UROMOL](https://github.com/sialindskrog/classifyNMIBC) classifiers. |
| **Exploratory_analysis.Rmd** | Runs quality control and PCA, and outputs cleaned counts/metadata tables used downstream. |
| **Sankey.R and Upset.R** | Create Sankey and UpsetPlots to visualize MIBC class and histological subtype distribution and co-occurences within patients. | 
| **Differential expression** | Identifies differentially expressed genes between groups using Limma + voom, with volcano plot, heatmap, and boxplot visualizations. |
| **GSEA** | Runs Gene Set Enrichment Analysis on DE genes against Hallmark and GO pathway sets. |
| **Microenvironment deconvolution** | Estimates immune/stromal cell composition per sample using MCPcounter via the [ImmuneDeconv interface](https://github.com/omnideconv/immunedeconv). |
| **Signature_scoring.Rmd** | Computes scores for given transcriptomic signatures. |


