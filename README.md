# Lewy Body Disease DNA Methylation Analysis

This repository contains the scripts, results, and full documentation associated with my Master's thesis titled:

**"Application of Supervised Machine Learning for Predictive Modeling of Lewy Body Disorders Using DNA Methylation Profiles"**  
by **Carlos C. Ureña Mateo** (University of Murcia, MSc Bioinformatics, 2025).

## Overview

This project implements a **reproducible and modular computational pipeline** for the analysis of DNA methylation data, aimed at identifying epigenetic patterns with discriminative value through the integration of prior biological knowledge and supervised machine learning methods.

The pipeline combines **clinically curated gene panels** with complementary supervised classification algorithms — **Support Vector Machines (SVM)** and **Random Forests (RF)** — to reduce dimensionality in a biologically informed manner and improve model interpretability and robustness.

As a case study, the pipeline is applied to DNA methylation data from patients affected by **Lewy body–related neurodegenerative diseases**, including **Parkinson’s disease (PD)**, **Parkinson’s disease dementia (PDD)**, and **dementia with Lewy bodies (DLB)**, together with healthy control samples. This design enables the exploration of disease-specific epigenetic signatures while accounting for shared molecular and clinical characteristics across related disorders.

For a detailed description of the methodology, analyses, and complete results, see the full project report in Spanish: [`LewyBody_Project_ES.pdf`](LewyBody_Project_ES.pdf) or in English: [`LewyBody_Project_EN.pdf`](LewyBody_Project_EN.pdf).

---

## Modules

Each module is implemented as an independent script located in the `scripts/` directory, allowing flexible reuse or selective execution depending on the analysis requirements.

### 1. `train_test_split.ipynb`
Splits normalized DNA methylation data and corresponding clinical metadata into **training (70%) and testing (30%) sets**, using stratification by diagnostic group and a fixed random seed for reproducibility.

**Input**
- Normalized methylation matrix (samples × probes)
- Clinical metadata table

**Output**
- Training methylation matrix
- Testing methylation matrix
- Corresponding clinical metadata files

---

### 2. `extract_gene_panels.R`
Extracts high-confidence disease-associated gene panels for targeted probe filtering.

**Input**
- PanelApp TSV files (Genomics England)
- Curated CSV gene panel files

**Output**
- Gene symbol vectors containing high-evidence genes

---

### 3. `filter_methylation_by_genes.R`
Filters DNA methylation probes based on overlap with clinically relevant gene panels.

**Input**
- Methylation matrix
- Probe annotation file
- Gene panel list

**Output**
- Filtered methylation matrix containing panel-associated probes

---

### 4. `prepare_for_ml_global.R`
Prepares datasets for binary machine learning classification at the global gene level.

**Input**
- Filtered methylation matrix
- Clinical metadata

**Output**
- ML-ready training and testing files per contrast

---

### 5. `train_models_all.R`
Trains and evaluates binary classification models using supervised machine learning.

**Implementation details**
- Implemented using the **caret** package
- Five-fold cross-validation
- Automatic hyperparameter tuning (`tuneLength = 5`)
- Algorithms:
  - Support Vector Machine (linear kernel)
  - Support Vector Machine (radial kernel)
  - Random Forest

**Evaluation metrics**
- Accuracy
- Cohen’s Kappa
- No Information Rate (NIR)
- McNemar’s test p-value
- Model p-value

**Input**
- Training datasets
- Model configuration parameters

**Output**
- Trained models
- Cross-validated performance metrics

---

### 6. `split_methylation_by_loci.R`
Splits filtered probes into biologically meaningful genomic loci: Promoter, TSS200, TSS1500, Body, 1stExon, 5’UTR, 3’UTR, and ExonBnd.

**Input**
- Filtered methylation matrix
- Probe genomic annotations

**Output**
- Separate methylation matrices per genomic locus

---

### 7. `prepare_for_ml_loci.R`
Prepares ML datasets independently for each genomic locus.

**Input**
- Locus-specific methylation matrices
- Clinical metadata

**Output**
- Training and testing datasets per locus and contrast

---

### 8. `train_models_loci.R`
Trains and evaluates binary classification models for each genomic locus and clinical contrast.

**Function**
- Fits independent models per locus to assess the predictive contribution of specific genomic regions.

**Implementation details**
- Implemented using the **caret** package
- Five-fold cross-validation
- Automatic hyperparameter tuning (`tuneLength = 5`)
- Same training strategy and evaluation framework as `train_models_all()`

**Algorithms**
- Support Vector Machine (linear kernel)
- Support Vector Machine (radial kernel)
- Random Forest

**Evaluation metrics**
- Accuracy
- Cohen’s Kappa
- No Information Rate (NIR)
- McNemar’s test p-value
- Model p-value

**Input**
- Locus-specific training datasets

**Output**
- Trained models
- Cross-validated performance metrics per locus and clinical contrast

---

## Results

All performance metrics and model outputs are stored in the `results/` directory, organized by **panel**, **locus**, **contrast**, and **model configuration**. The directory contains structured PDF reports summarizing:

- Accuracy, Balanced Accuracy, Cohen’s Kappa, ROC AUC, and p-values
- Locus-specific and global model performance
- Covariate-adjusted Random Forest analyses

---

## Key Findings

- **Global models**: No statistically significant results after multiple testing correction.  
- **Locus-specific analyses**: Significant signals detected in promoter and TSS200 regions of the Adult Neurodegenerative Disorders panel, particularly for **PD vs CTRL** using radial SVM.  
- **Covariate-adjusted Random Forest models**: Significant signals only in negative control panels, suggesting demographic imbalances may affect associations.  
- **Conclusion**: Incorporating biological prior knowledge and functional genomic stratification is critical for modeling high-dimensional DNA methylation data effectively.

---

## License

This repository is provided for academic and research purposes. Please cite the thesis if used in your work.


