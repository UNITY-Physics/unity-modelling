# UNITY Network Modeling Repository

Welcome to the **UNITY Network Modeling Repository** 👋  

This repository is part of the [UNITY Network](https://www.unity-mri.com/), a global initiative to improve child health through the use of **ultra-low-field MRI** and other accessible neuroimaging tools. Here, researchers and developers collaborate to share, test, and improve **code for modeling data** across diverse populations and contexts.

---

## 🌍 Purpose

The UNITY Network spans many research sites worldwide, working together to ensure that advances in neuroimaging are **globally inclusive** and **locally relevant**.  

This repository provides a **shared workspace** where contributors can:  
- Upload and maintain code for data modeling cognitive outcomes.  
- Share best practices and examples for reproducible research.  
- Collaborate on common tools and pipelines, reducing duplication of effort across sites.  

---

## 📂 Repository Structure

```
unity-modeling/
│
├── NeuroCogRisk.R              # Main analysis script for neurodevelopmental risk prediction
├── models/                    # Modeling functions and data preparation scripts
│   ├── xgboost.R              # XGBoost modeling functions with SHAP analysis
│   └── data_prep/             # Data preparation scripts for different cohorts
│       ├── flywheel_eeg_merge.R           # Merge Khula and SPACE EEG data to FlyWheel MRI dataset
│       ├── flywheel_khula_eeg.R          # Create Khula EEG dataset for FlyWheel
│       ├── flywheel_space_eeg.R           # Create SPACE EEG dataset for FlyWheel
│       ├── flywheel_terminal_gsed_bsid.R # Extract latest outcome scores per subject
│       └── superfield_eeg_other_var_data_prep.R # Prepare MRI, EEG, and demographic data for analysis
├── utils/                     # Utility functions
│   ├── preprocessing.R        # Main data preprocessing pipeline for GATES cohort
│   └── VennDiagram.R          # Create Venn diagrams for cross-cohort comparisons
├── data/                      # Processed datasets (output from preprocessing)
├── notebooks/                 # Jupyter notebooks for exploration, tutorials, or demos
├── docs/                      # Documentation and guidelines
└── tests/                     # Unit and integration tests
```

---

## 🔬 Main Components

### Analysis Scripts

#### `NeuroCogRisk.R`
Main analysis script that fits XGBoost models to predict neurodevelopmental outcomes (BSID scores and DAZ/GSED scores) using early-life predictors. The script:
- Analyzes multiple cohorts: PRIMES, Khula+SPACE, MINE, PRISMA-Kenya, PRISMA-Zambia
- Tests different predictor combinations: baseline demographics, brain volumes (main field and spinal cord), and EEG features
- Generates SHAP plots and performance summaries for model interpretability
- Saves results as PDF files with performance tables

#### `models/xgboost.R`
Core modeling functions including:
- `fit_xgboost_miss()`: Fits XGBoost models with automatic parameter tuning and missing data imputation
- `tune_xgboost_params()`: Automatically selects optimal parameters to minimize overfitting
- `create_shap_page_xgb()`: Creates SHAP (SHapley Additive exPlanations) visualizations
- `perf_row()`: Generates performance summary tables

### Data Preparation Scripts

#### `utils/preprocessing.R`
Main preprocessing pipeline that:
- Cleans and standardizes GATES cohort data
- Extracts BSID (Bayley Scales) and DAZ scores
- Calculates brain volume residuals (head circumference-adjusted)
- Extracts EEG frequency band features
- Creates analysis-ready datasets for specific cohorts

#### `models/data_prep/`
Cohort-specific data preparation scripts:
- **flywheel_eeg_merge.R**: Merges EEG data from Khula and SPACE cohorts with MRI data from FlyWheel
- **flywheel_khula_eeg.R**: Processes Khula cohort EEG data (spectral power, microstates, total power)
- **flywheel_space_eeg.R**: Processes SPACE cohort EEG data
- **flywheel_terminal_gsed_bsid.R**: Extracts the latest available outcome scores (GSED/BSID) per subject
- **superfield_eeg_other_var_data_prep.R**: Prepares combined datasets with MRI, EEG, and demographic variables for predictive analysis

### Utility Scripts

#### `utils/VennDiagram.R`
Creates Venn diagrams to visualize:
- Overlap of brain regions measured across different cohorts (UCT, MINE, PRISMA, PRIMES)
- Overlap of baseline variables available across studies