# MSc Bioinformatics Thesis – Lewy Body Disease DNA Methylation Analysis

This repository contains the scripts and results associated with my Master's thesis titled:

**"Application of Supervised Machine Learning for Predictive Modeling of Lewy Body Disorders Using DNA Methylation Profiles"**  
by **Carlos C. Ureña Mateo** (University of Murcia, MSc Bioinformatics, 2025).

## Overview

This project implements a **reproducible and modular computational pipeline** for the analysis of DNA methylation data, aimed at identifying epigenetic patterns with discriminative value through the integration of prior biological knowledge and supervised machine learning methods.

The pipeline combines **clinically curated gene panels** with complementary supervised classification algorithms — **Support Vector Machines (SVM)** and **Random Forests (RF)** — to reduce dimensionality in a biologically informed manner and improve model interpretability and robustness.

As a case study, the pipeline is applied to DNA methylation data from patients affected by **Lewy body–related neurodegenerative diseases**, including **Parkinson’s disease (PD)**, **Parkinson’s disease dementia (PDD)**, and **dementia with Lewy bodies (DLB)**, together with healthy control samples. This design enables the exploration of disease-specific epigenetic signatures while accounting for shared molecular and clinical characteristics across related disorders.

## Modules

Each module is implemented as an independent script located in the `scripts/` directory, allowing flexible reuse or selective execution depending on the analysis requirements.

### 1. `train_test_split.ipynb`

Splits normalized DNA methylation data and corresponding clinical metadata into independent training and testing sets.

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

### 5. `train_models_all()`

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

Splits filtered probes into biologically meaningful genomic loci: Promoter, TSS200, TSS1500, Body, 1stExon, 5’UTR, 3’UTR and ExonBnd.

**Input**
- Filtered methylation matrix
- Probe genomic annotations

**Output**
- Separate methylation matrices per genomic locus

---

### 7. `prepare_for_ml_loci()`

Prepares ML datasets independently for each genomic locus.

**Input**
- Locus-specific methylation matrices
- Clinical metadata

**Output**
- Training and testing datasets per locus and contrast

---
### 8. `train_models_loci()`

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

## Methodology (Case od study)

This project implements a **modular, reproducible DNA methylation analysis pipeline**, designed to identify epigenetic patterns with discriminative value through integration of prior biological knowledge and supervised machine learning.

**Data preparation**  
Normalized methylation data and clinical metadata were split into **training (70%) and testing (30%) sets**, using stratification by diagnostic group and a fixed random seed (`train_test_split.ipynb`). Clinically relevant gene panels were obtained from **Genomics England PanelApp** and curated CSV sources, retaining only **high-evidence genes** (`get_green_genes_from_genomics_england()`, `get_high_evidence_genes_from_csv()`).

**Panel-based probe filtering**  
CpG probes overlapping panel genes were extracted to reduce dimensionality and focus on disease-relevant regions (`filter_methylation_by_genes()`).

**Global panel-level modeling**  
Binary classification models were trained for each panel and contrast (PD vs CTRL, PDD vs CTRL, DLB vs CTRL, neuro vs CTRL) using all retained probes (`prepare_for_ml_global()`, `train_models_all()`). Algorithms: **SVM (linear and radial)** and **Random Forest**, with 5-fold stratified cross-validation and automatic hyperparameter tuning. Performance was evaluated on independent test sets.

**Locus-specific modeling**  
Probes were grouped into **functional genomic regions** (promoter, TSS200, TSS1500, body, UTRs, exon boundaries) (`split_methylation_by_loci()`). Models were trained per **panel × locus × contrast** (`prepare_for_ml_loci()`, `train_models_loci()`), enabling detection of localized epigenetic signals.

**Adjustment for clinical covariates**  
Age and sex were incorporated into locus-level Random Forest models to assess their effect on model performance.

**Results**  
Performance metrics (Accuracy, Balanced Accuracy, Kappa, ROC AUC, p-values) are stored in a **structured results directory** organized by panel, locus, contrast, and model configuration.

---

## Key Findings

- Global models did not yield statistically significant results after multiple testing correction.  
- Locus-specific analyses revealed significant signals in promoter and TSS200 regions of the Adult Neurodegenerative Disorders panel, particularly for Parkinson’s disease vs Control (PD vs CTRL) using radial SVM.  
- Random Forest models including covariates detected significant signals only in negative control panels, suggesting that demographic imbalances may drive associations rather than disease-specific methylation differences.  
- These findings highlight the importance of incorporating biological prior knowledge and functional genomic stratification when modeling high-dimensional DNA methylation data.


