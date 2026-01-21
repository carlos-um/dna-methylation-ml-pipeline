# Lewy Body Disease DNA Methylation Analysis

## Overview

This project implements a **reproducible and modular computational pipeline** for the analysis of DNA methylation data, aimed at identifying epigenetic patterns with discriminative value through the integration of prior biological knowledge and supervised machine learning methods.

The pipeline combines **clinically curated gene panels** with complementary supervised classification algorithms — **Support Vector Machines (SVM)** and **Random Forests (RF)** — to reduce dimensionality in a biologically informed manner and improve model interpretability and robustness.

As a case study, the pipeline is applied to DNA methylation data from patients affected by **Lewy body–related neurodegenerative diseases**, including **Parkinson’s disease (PD)**, **Parkinson’s disease dementia (PDD)**, and **dementia with Lewy bodies (DLB)**, together with healthy control samples. This design enables the exploration of disease-specific epigenetic signatures while accounting for shared molecular and clinical characteristics across related disorders.

For detailed methodology, comprehensive analyses, and complete results, see:
- Full thesis report: [`report/LewyBody_Project_ES.pdf`](report/LewyBody_Project_ES.pdf)
- Pipeline implementation: [`scripts/`](scripts/)
- Analysis results: [`results/`](results/)


## Key Findings

- **Global models**: No statistically significant results after multiple testing correction.  
- **Locus-specific analyses**: Significant signals detected in promoter and TSS200 regions of the Adult Neurodegenerative Disorders panel, particularly for **PD vs CTRL** using radial SVM.  
- **Covariate-adjusted Random Forest models**: Significant signals only in negative control panels, suggesting demographic imbalances may affect associations.  
- **Conclusion**: Incorporating biological prior knowledge and functional genomic stratification is critical for modeling high-dimensional DNA methylation data effectively.

---

## Author

**Carlos C. Ureña Mateo**  
MSc Bioinformatics – University of Murcia

---

## License

This repository is provided for academic and research purposes. Please cite the thesis if used in your work.


