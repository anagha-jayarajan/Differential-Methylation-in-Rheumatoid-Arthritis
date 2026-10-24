
# DNA Methylation Analysis in Rheumatoid Arthritis

## Project Overview

This project investigates genome-wide DNA methylation differences between rheumatoid arthritis (RA) patients and healthy controls using publicly available DNA methylation data from the NCBI Gene Expression Omnibus (GEO) database.

The analysis uses the **GSE42861 dataset**, generated using the Illumina HumanMethylation450 BeadChip platform, to identify differentially methylated positions (DMPs), explore disease-associated methylation patterns, evaluate sample quality, and investigate potential epigenetic changes associated with rheumatoid arthritis.

Three complementary R-based bioinformatics pipelines—**ChAMP, minfi, and SeSAMe**—were used for data preprocessing, quality control, normalization, and methylation analysis.

## Objectives

- Analyze genome-wide DNA methylation patterns in rheumatoid arthritis.
- Assess the quality and reliability of raw IDAT intensity data.
- Perform methylation preprocessing and normalization using established Bioconductor tools.
- Identify differentially methylated CpG positions between RA cases and healthy controls.
- Visualize disease-associated methylation patterns using PCA, MDS, volcano plots, and hierarchical clustering heatmaps.
- Investigate the genomic distribution and potential biological relevance of differentially methylated genes.
- Compare complementary methylation analysis workflows and document their limitations.

## Dataset Information

| Parameter | Details |
|---|---|
| GEO accession | GSE42861 |
| Organism | *Homo sapiens* |
| Disease | Rheumatoid arthritis |
| Tissue | Peripheral blood leukocytes (PBLs) |
| Platform | Illumina HumanMethylation450 BeadChip |
| GEO platform accession | GPL13534 |
| Samples analyzed | 10 |
| RA cases | 5 |
| Healthy controls | 5 |
| Data format | Illumina IDAT files |
| Analysis language | R |

**Dataset:** [GSE42861 – NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42861)

The analysis uses a 10-sample subset of the larger GEO study. Consequently, the findings should be interpreted as exploratory rather than as definitive evidence of rheumatoid arthritis-associated epigenetic changes.

## Tools and Technologies

### R and Bioconductor

- **R:** Statistical computing and data analysis.
- **ChAMP:** Methylation preprocessing, quality control, normalization, and differential methylation analysis.
- **minfi:** IDAT processing, sample-level quality control, normalization, and differential methylation analysis.
- **SeSAMe:** IDAT-level processing, signal correction, and methylation quality assessment.
- **ggplot2:** Data visualization.
- **pheatmap:** Hierarchical clustering and methylation heatmaps.

## Analysis Workflow

### 1. Raw Data Processing

The raw Illumina methylation data were obtained in IDAT format. Each sample was represented by red and green channel intensity files.

Raw intensity summaries were examined to evaluate signal quality and identify potential technical variation.

### 2. Quality Control

Quality control was performed to assess the reliability of the methylation measurements.

The workflow included:

- Detection p-value assessment.
- Methylated and unmethylated signal intensity evaluation.
- Illumina control probe assessment.
- Sample-level quality assessment using minfi.
- Examination of raw beta-value density distributions.
- Identification of possible array-related batch effects.

All ten samples passed the reported mean detection p-value quality threshold.

### 3. ChAMP Pipeline

ChAMP was used as the primary differential methylation analysis pipeline.

Key processing steps included:

- Loading IDAT files and generating methylation beta values.
- Filtering probes with detection p-values greater than 0.01.
- Removing SNP-associated, multi-hit, and non-CpG probes.
- BMIQ (Beta Mixture Quantile) normalization.
- Singular Value Decomposition (SVD) to investigate technical and biological variation.
- Differentially Methylated Position (DMP) analysis using empirical Bayes moderation.

### 4. minfi Pipeline

The minfi workflow provided complementary quality control and preprocessing.

The analysis included:

- Importing IDAT files using `read.metharray.exp()`.
- Assessing sample quality using `plotQC()`.
- Evaluating Illumina control probe performance.
- Applying `preprocessQuantile()` normalization.
- Examining methylation density and multidimensional scaling plots.
- Performing differential methylation analysis using `dmpFinder()`.

### 5. SeSAMe Pipeline

SeSAMe was applied to the raw IDAT data using `openSesame()` with p-value masking.

The pipeline was used primarily to assess signal quality and provide a complementary view of the methylation data. Its differential methylation output contained NA p-values, preventing independent statistical validation of the DMP findings in this subset.

### 6. Exploratory Data Analysis and Visualization

The workflow included several visualizations to understand sample-level patterns and disease-associated methylation differences:

- Raw beta-value density plots.
- minfi methylated versus unmethylated intensity plots.
- Multidimensional scaling (MDS).
- Singular Value Decomposition (SVD) analysis.
- Principal Component Analysis (PCA).
- Volcano plots of differential methylation results.
- Hierarchical clustering heatmaps of the top 50 DMPs.

These visualizations were used to examine sample clustering, identify potential batch effects, and explore differences between RA cases and controls.

## Key Results

### Differential Methylation

The ChAMP analysis reported the following results:

| Metric | Reported result |
|---|---:|
| Total probes tested | 3,044 |
| Significant DMPs (adjusted p-value < 0.05) | 3,044 |
| Hypomethylated DMPs (delta beta < -0.1) | 2,053 |
| Hypermethylated DMPs (delta beta > 0.1) | 560 |
| Top reported DMP | cg21193926 |
| Associated gene | TGFB3 |
| Top reported adjusted p-value | 6.73 × 10⁻⁷ |
| Maximum absolute delta beta | 0.87 |
| Probe associated with maximum absolute delta beta | cg27286614 (MAD1L1) |

The reported results suggest widespread methylation differences between the two groups. However, the reported DMP counts and effect-direction summaries require further validation before being interpreted as a final set of statistically reliable findings.

### Notable Genes

| Gene | Reported methylation pattern | Biological context |
|---|---|---|
| TGFB3 | Hypermethylation | TGF-β signaling and immune regulation |
| MAD1L1 | Hypomethylation | Mitotic checkpoint and cell-cycle regulation |
| BRDT | Hypermethylation | Chromatin regulation |
| TLE1 | Hypomethylation | Wnt signaling and transcriptional regulation |
| OTUD5 | Hypomethylation | Deubiquitination and NF-κB-related signaling |
| TSPY1 | Hypermethylation | Y-chromosome-associated biology |

These genes represent candidates for further investigation. Their methylation patterns alone do not establish changes in gene expression or prove a causal role in rheumatoid arthritis.

### Genomic Distribution

The reported differentially methylated positions were annotated across several genomic regions, including:

- Gene bodies.
- Transcription start site regions (TSS200 and TSS1500).
- First exons.
- 5′ and 3′ untranslated regions.
- Intergenic regions.
- CpG islands, shores, shelves, and open-sea regions.

The report describes enrichment of CpG-island-associated DMPs, suggesting that regulatory genomic regions may be affected. This observation requires confirmation using validated probe annotations and independently reproducible statistical results.

### Sample Clustering and Batch Effects

PCA and MDS were used to explore the relationship between methylation patterns and disease status.

The report describes partial or distinct separation between RA and control samples, alongside potential technical variation associated with array slides. In particular, one control sample originated from a different Sentrix array.

These observations highlight the importance of evaluating batch effects and sample-level confounding in methylation studies.

## Limitations

- **Small sample size:** The analysis includes only five RA cases and five controls, limiting statistical power and generalizability.
- **Subset analysis:** Results from ten samples may not represent the full GSE42861 cohort.
- **Batch effects:** Differences between array slides may influence observed methylation patterns.
- **SeSAMe statistical limitations:** NA p-values prevented independent DMP validation using this pipeline.
- **Result validation:** The reported number of significant DMPs and the methylation effect-direction summaries require rechecking against the underlying analysis outputs.
- **Biological interpretation:** Differential methylation does not directly demonstrate altered gene expression, functional consequences, or causality.

The findings should therefore be considered preliminary until the complete preprocessing, statistical testing, and validation workflow has been independently reproduced.

## Reproducibility and Future Work

Potential next steps include:

1. Reanalyzing the complete available GSE42861 cohort.
2. Verifying probe filtering, normalization, and sample metadata.
3. Accounting for batch effects and relevant biological covariates.
4. Rechecking adjusted p-values, delta beta values, and methylation direction.
5. Comparing DMP results across ChAMP, minfi, and SeSAMe where statistically supported.
6. Performing gene-level and pathway enrichment analyses.
7. Integrating methylation findings with gene-expression data where suitable data are available.
8. Validating candidate CpG sites and genes in independent datasets.

## Repository Structure

A suggested repository structure is:

```text
rheumatoid-arthritis-dna-methylation-analysis/
├── README.md
├── data/
│   └── README.md
├── scripts/
│   ├── champ_analysis.R
│   ├── minfi_analysis.R
│   └── sesame_analysis.R
├── results/
│   ├── qc/
│   ├── differential_methylation/
│   └── figures/
├── report/
│   └── DNA_Methylation_Report1.pdf
└── .gitignore
```

This is a suggested organization; include only files and scripts that actually exist in the repository.

## Software Requirements

- R version 4.0 or later, as reported in the project documentation.
- Bioconductor packages: ChAMP, minfi, and SeSAMe.
- Visualization packages: ggplot2 and pheatmap.

Package versions and installation instructions should be documented to support reproducibility.

## References

1. NCBI Gene Expression Omnibus. **GSE42861: Differential DNA methylation in Rheumatoid Arthritis.**  
   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42861

2. ChAMP: Chip Analysis Methylation Pipeline.  
   https://bioconductor.org/packages/ChAMP/

3. minfi: Analyze Illumina Infinium DNA methylation arrays.  
   https://bioconductor.org/packages/minfi/

4. SeSAMe: SEnsitive Step-wise Analysis of Methylation data.  
   https://bioconductor.org/packages/sesame/

## Authors

- Anagha Jayarajan


## Acknowledgements

The authors acknowledge the NCBI Gene Expression Omnibus for providing publicly accessible DNA methylation data and the Bioconductor community for developing the computational tools used in this project.

provider's terms.

