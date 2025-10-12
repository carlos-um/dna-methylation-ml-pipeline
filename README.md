# MSc Bioinformatics Thesis – Lewy Body Disease DNA Methylation Analysis

This repository contains the scripts, and results associated with my Master's thesis titled:

**"Application of Supervised Machine Learning for Predictive Modeling of Lewy Body Disorders Using DNA Methylation Profiles"**  
by **Carlos C. Ureña Mateo** (University of Murcia, MSc Bioinformatics, 2025).

## Overview

Lewy body diseases, including Parkinson’s disease, Parkinson’s disease dementia, and dementia with Lewy bodies, are neurodegenerative disorders that may share epigenetic alterations, particularly in DNA methylation profiles. The aim of this thesis was to design and implement a reproducible and modular computational pipeline to analyze DNA methylation data, integrating clinically curated gene panels with supervised machine learning models to identify biologically meaningful and discriminative epigenetic patterns.

## Methods

- **Data preprocessing:** Filtering methylation arrays using 16 biologically-informed gene panels (including four related to neurodegeneration and several control panels).  
- **Machine learning:** Classification models were trained using SVM (radial and linear kernels) and Random Forest (ranger), both at a global level (entire panel) and locus-specific (e.g. promoters, TSS regions).  
- **Covariate analysis:** Additional models included clinical covariates (age and sex) to assess potential confounding effects.  
- **Pipeline design:** Fully modular, reproducible workflow implemented in R and Python, compatible with HPC environments.

## Repository Structure

- **`scripts/`** – All R and Python scripts used to process data, perform filtering, and train models.  
- **`results/`** – Output files and results from model training and statistical analyses.  

## Key Findings

- Global models did not yield statistically significant results after multiple testing correction.  
- Locus-specific analysis revealed significant signals in promoter and TSS200 regions of the Adult Neurodegenerative Disorders panel, especially for Parkinson’s disease vs Control (PD vs CTRL) using radial SVM.  
- Random Forest models including covariates detected significant signals only in negative control panels, suggesting demographic imbalances may drive associations rather than disease-specific methylation differences.  
- These findings highlight the importance of integrating biological knowledge and functional genomic stratification in methylation studies.

## Contact

**Carlos C. Ureña Mateo**  
Email: cc.urenamateo@gmail.com

