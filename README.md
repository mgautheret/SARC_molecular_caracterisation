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

Overall, the pipeline takes as input raw RNA-seq counts, TPM counts and sample metadata, and produces molecular classifying, cleaned tables, differential expression results, and various QC/analysis plots.

| Script | Description |
|---|---|
| **Molecular_classification.Rmd** | Classifies samples into MIBC/NMIBC subtypes using the [consensusMIBC](https://github.com/cit-bioinfo/consensusMIBC) and [UROMOL](https://github.com/sialindskrog/classifyNMIBC) classifiers. |
| **QC_Exploratory_analysis.Rmd** | Runs quality control and PCA, and outputs cleaned counts/metadata tables used downstream. |
| **Sankey.R and Upset.R** | Create Sankey and UpsetPlots to visualize MIBC class and histological subtype distribution and co-occurences within patients. | 
| **Differential_expression_limma** | Identifies differentially expressed genes between groups using Limma + voom, with volcano plot, heatmap, and boxplot visualizations. |
| **GSEA.Rmd** | Runs Gene Set Enrichment Analysis on DE genes against Hallmark and GO pathway sets. |
| **ImmuneDeconv.Rmd** | Estimates immune/stromal cell composition per sample using MCPcounter via the [ImmuneDeconv interface](https://github.com/omnideconv/immunedeconv). |
| **Signature_scoring.Rmd** | Computes scores for given transcriptomic signatures. |

## Dependencies 
- readxl (1.5.0)
- UCell (2.14.0)
- GOSemSim (2.36.0)
- org.Hs.eg.db (3.22.0)
- msigdbdf (25.1.0)
- msigdbr (26.1.0)
- enrichplot (1.30.5)
- clusterProfiler (4.18.4)
- classifyNMIBC (1.1.0)
- consensusMIBC (1.1.0)
- immunedeconv (2.1.4)
- ComplexHeatmap (2.26.1)
- edgeR (4.8.2)
- networkD3 (0.4.1)
- UpSetR (1.4.1)

## Acknowledgment 
Thanks to Guerric Gilbert (Institut Curie IMPACT team, Paris) for bioinformatics support.
- Kamoun, A et. al. A Consensus Molecular Classification of Muscle-invasive Bladder Cancer, Eur Urol (2019), doi: https://doi.org/10.1016/j.eururo.2019.09.006
- Lindskrog and Prip et al. An integrated multi-omics analysis identifies prognostic molecular subtypes of non-muscle-invasive bladder cancer. Nat Commun. 2021. PMID: 33863885. DOI: 10.1038/s41467-021-22465-w
- Sturm G, Finotello F, Petitprez F, Zhang JD, Baumbach J, Fridman WH, List M, Aneichyk T. Comprehensive evaluation of transcriptome-based cell-type quantification methods for immuno-oncology. Bioinformatics. 2019 Jul 15;35(14):i436-i445. doi: 10.1093/bioinformatics/btz363.

## Author 
Madeleine Gautheret, Master's student and intern, Université de Rennes/IMPACT team


  
