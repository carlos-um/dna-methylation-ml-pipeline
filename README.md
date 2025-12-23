# MSc Bioinformatics Thesis – Lewy Body Disease DNA Methylation Analysis

This repository contains the scripts, and results associated with my Master's thesis titled:

**"Application of Supervised Machine Learning for Predictive Modeling of Lewy Body Disorders Using DNA Methylation Profiles"**  
by **Carlos C. Ureña Mateo** (University of Murcia, MSc Bioinformatics, 2025).

## Overview

Lewy body diseases — including Parkinson’s disease (PD), Parkinson’s disease dementia (PDD), and dementia with Lewy bodies (DLB) — are neurodegenerative disorders for which disease-specific epigenetic alterations remain incompletely characterized. DNA methylation provides a high-resolution molecular layer to investigate whether these conditions can be distinguished from neurologically healthy controls based on epigenetic profiles.

In this thesis, a reproducible and modular bioinformatics pipeline was developed to analyze DNA methylation data using a hypothesis-driven strategy. Clinically curated gene panels were incorporated at the initial stages of the pipeline as prior biological knowledge to guide probe selection and reduce dimensionality before any statistical modeling.

Following this biologically informed filtering step, supervised machine learning models were applied to evaluate the discriminative capacity of panel-specific methylation profiles. Models were trained both at a global panel level and after stratification by functional genomic regions, enabling the detection of localized epigenetic signals associated with disease status rather than assuming shared methylation patterns across disorders.

## Methods

- **Data preprocessing:** Biologically informed probe filtering of preprocessed DNA methylation data using 16 curated gene panels (including four related to neurodegeneration and several control panels).  
- **Machine learning:** Classification models were trained using SVM (radial and linear kernels) and Random Forest (ranger), both at a global panel level and in locus-specific analyses (e.g. promoters, TSS regions).  
- **Covariate analysis:** Additional models included clinical covariates (age and sex) to assess potential confounding effects.  
- **Pipeline design:** Fully modular, reproducible workflow implemented in R and Python, compatible with HPC environments.

## Key Findings

- Global models did not yield statistically significant results after multiple testing correction.  
- Locus-specific analyses revealed significant signals in promoter and TSS200 regions of the Adult Neurodegenerative Disorders panel, particularly for Parkinson’s disease vs Control (PD vs CTRL) using radial SVM.  
- Random Forest models including covariates detected significant signals only in negative control panels, suggesting that demographic imbalances may drive associations rather than disease-specific methylation differences.  
- These findings highlight the importance of incorporating biological prior knowledge and functional genomic stratification when modeling high-dimensional DNA methylation data.


