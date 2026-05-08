# ============================================================================
# GATES Data Preprocessing Script
# Date: September 25, 2025
# 
# Purpose: This script preprocesses GATES cohort data for analysis, including:
#   - Data cleaning and variable transformation
#   - Extraction of BSID (Bayley Scales) and DAZ scores
#   - Brain volume measurements and residual calculations
#   - EEG feature extraction
#   - Creation of analysis-ready datasets for specific cohorts
# ============================================================================

# Load required packages
library(dplyr)      # For data manipulation and piping
library(tidyr)      # For data tidying functions
library(purrr)      # For functional programming tools
library(ggplot2)    # For plotting (if needed)
library(patchwork)  # For combining plots (if needed)

# ============================================================================
# DATA IMPORT AND INITIAL CLEANING
# ============================================================================

# Read the main GATES dataset
# This contains demographic, anthropometric, and developmental assessment data
# Note: Update the file path to match your local data location
threshold <- 1  # Age threshold in years for early predictor selection
gates <- read.csv('/Users/yidong/Google Drive/Projects/Data/GATES/UNITY-dataTemplate_MRIEEG_clean.csv')

# Standardize PRISMA cohort naming by appending country location
# This ensures PRISMA cohorts from different countries are distinguishable
gates <- gates %>%
  mutate(CohortName = if_else(
    CohortName == "PRISMA",
    paste(CohortName, CohortLocation_country, sep = "-"),
    CohortName
  ))

# Clean the gates dataset in place: convert empty strings to NA for key character variables
# Also perform data transformations and create derived variables
gates <- gates %>%
  mutate(
    # Basic identifiers and demographics
    id = UniqueStudyID,
    cohort = CohortName,
    sex = ifelse(childBiologicalSex == "", NA, childBiologicalSex),
    
    # Age calculation: use months if days is missing, convert to days
    age = ifelse(is.na(childTimepointAge_days) & !is.na(childTimepointAge_months), 
                 childTimepointAge_months * 365.25 / 12, 
                 childTimepointAge_days),
    age_eeg = eeg_age_days,
    
    # Developmental scores
    daz = MAX_gsed_LongForm_DAZScore,
    age_daz = MAX_gsed_LongForm_DAZScore_Age_days,
    
    # BSID scores (Bayley Scales of Infant Development)
    bsidc = MAX_bsid_CCS,      # Cognitive composite score
    age_bsidc = MAX_bsid_CCS_Age_days,
    bsidl = MAX_bsid_LCS,      # Language composite score
    age_bsidl = MAX_bsid_LCS_Age_days,
    bsidm = MAX_bsid_MCS,      # Motor composite score
    age_bsidm = MAX_bsid_MCS_Age_days,
    
    # Maternal and birth characteristics
    medu = MaternalEducation_schoolingYears,           # Maternal education in years
    ga = childGestation_weeks,                         # Gestational age in weeks
    ma = maternalAgeAtChildBirth_years,                # Maternal age at childbirth
    fs = timepointFamilySize,
    
    # Anthropometric measurements
    hcmri = childTimepointHC_MRI_cm,                   # Head circumference from MRI
    hc = childTimepointHC_measured_cm,
    bw = childBirthWeight_kgs,                         # Birth weight in kg
    bl = childBirthLength_inches * 2.54,               # Birth length converted to cm
    
    # Growth indicators (Z-scores for anthropometric measurements)
    waz = WAZ,                                          # Weight-for-age Z-score
    haz = HAZ,                                          # Height-for-age Z-score
    whz = WHZ,                                          # Weight-for-height Z-score
    
    # Brain volume measurements (in cm³)
    # Note: Original values are in mm³, dividing by 10000 converts to cm³
    # Bilateral structures are averaged across left and right hemispheres
    caudate = (MM_left_caudate + MM_right_caudate) / 10000,
    putamen = (MM_left_putamen + MM_right_putamen) / 10000,
    thalamus = (MM_left_thalamus + MM_right_thalamus) / 10000,
    ventricles = MM_ventricles / 10000,
    anterior_callosum = MM_anterior_callosum / 10000,
    posterior_callosum = MM_posterior_callosum / 10000,
    central_callosum = (MM_mid_posterior_callosum + MM_central_callosum + MM_mid_anterior_callosum) / 10000,
    corpus_callosum = anterior_callosum + posterior_callosum + central_callosum,
    supratentorial_csf = MM_supratentorial_csf / 10000,
    supratentorial_tissue = MM_supratentorial_tissue / 10000,
    cerebellum = MM_cerebellum / 10000,
    globus_pallidus = (MM_left_globus_pallidus + MM_right_globus_pallidus) / 10000,
    
    # Spinal cord (SC) measurements (in cm³)
    # These are from a different imaging field/sequence
    sc_csf = sc_csf / 10000, 
    sc_gray_matter = sc_gray_matter / 10000, 
    sc_white_matter = sc_white_matter / 10000,
    sc_caudate = (sc_left_caudate + sc_right_caudate) / 10000,
    sc_lentiform = (sc_left_lentiform + sc_right_lentiform) / 10000,
    sc_hippocampus = (sc_left_hippocampus + sc_right_hippocampus) / 10000,
    sc_thalamus = (sc_left_thalamus + sc_right_thalamus) / 10000,
    sc_corpus_callosum = sc_corpus_callosum / 10000,
    
    # EEG frequency band power (averaged across channels)
    # These are mean values across all channels for each frequency band
    LowAlpha = rowMeans(select(., all_of(names(gates)[grepl("LowAlpha", names(gates))])), na.rm = TRUE),
    HighAlpha = rowMeans(select(., all_of(names(gates)[grepl("HighAlpha", names(gates))])), na.rm = TRUE),
    Beta = rowMeans(select(., all_of(names(gates)[grepl("Beta", names(gates))])), na.rm = TRUE),
    Gamma = rowMeans(select(., all_of(names(gates)[grepl("Gamma", names(gates))])), na.rm = TRUE), 
    # EEG power ratios (commonly used in neurodevelopmental research)
    ThetaBeta = rowMeans(select(., all_of(names(gates)[grepl("theta_beta_ratio", names(gates))])), na.rm = TRUE), 
    HighAlphaDelta = rowMeans(select(., all_of(names(gates)[grepl("high_alpha_delta_ratio", names(gates))])), na.rm = TRUE), 
    AlphaPeak = alpha_peak_binary  # Binary indicator for presence of alpha peak
  ) %>%
  # Select only the cleaned and derived variables
  select(id, cohort, sex, age, age_eeg, daz, bsidc, bsidl, bsidm, age_daz, age_bsidc, age_bsidl, age_bsidm,
         medu, ga, fs, ma, hcmri, hc, 
         bw, bl, waz, haz, whz,
         caudate, putamen, thalamus, ventricles,
         anterior_callosum, posterior_callosum, central_callosum, corpus_callosum,
         supratentorial_csf, supratentorial_tissue, cerebellum,
         globus_pallidus,
         sc_csf, sc_gray_matter, sc_white_matter,
         sc_caudate, sc_lentiform, sc_hippocampus, sc_thalamus, sc_corpus_callosum,
         LowAlpha, HighAlpha, Beta, Gamma, ThetaBeta, HighAlphaDelta, AlphaPeak) %>%
  # Remove rows with missing critical information
  drop_na(id)

# Fill missing baseline variables using repeated measurements from same child
# This imputes missing values for variables that shouldn't change over time
gates <- gates %>%
  group_by(id) %>%
  mutate(across(all_of(c("cohort", "sex", "medu", "ga", "fs", "ma", "bw", "bl")), 
                ~ifelse(is.na(.), first(na.omit(.)), .))) %>%
  ungroup() %>%
  drop_na(cohort, sex)

# Quality check: Identify children with inconsistent sex values across timepoints
# This is a data quality check - should return empty if data is consistent
# If non-empty, indicates data entry errors that need investigation
gates %>%
  group_by(id) %>%
  filter(n() > 1) %>%  # Only children with multiple rows
  filter(n_distinct(sex) > 1) %>%  # Different sex values across rows
  ungroup() %>%
  arrange(id, age) %>% 
  pull(id)

# ============================================================================
# EXTRACT BSID SCORES (Bayley Scales of Infant Development)
# ============================================================================
# Extract BSID scores separately for each domain (Cognitive, Language, Motor)
# Each dataset contains one score per child per assessment timepoint
# Multiple timepoints per child are preserved for longitudinal analysis

# BSID Cognitive Composite Score
bsidc <- gates %>%
  mutate(bsid = bsidc) %>%
  select(cohort, id, sex, age = age_bsidc, bsid) %>%
  drop_na(bsid) %>%  # Remove rows with missing BSID scores
  arrange(cohort, id, age) %>%
  distinct(id, age, .keep_all = TRUE)  # Remove duplicate assessments at same age

# BSID Language Composite Score
bsidl <- gates %>%
  mutate(bsid = bsidl) %>%
  select(cohort, id, sex, age = age_bsidl, bsid) %>%
  drop_na(bsid) %>%
  arrange(cohort, id, age) %>%
  distinct(id, age, .keep_all = TRUE)

# BSID Motor Composite Score
bsidm <- gates %>%
  mutate(bsid = bsidm) %>%
  select(cohort, id, sex, age = age_bsidm, bsid) %>%
  drop_na(bsid) %>%
  arrange(cohort, id, age) %>%
  distinct(id, age, .keep_all = TRUE)

# ============================================================================
# EXTRACT DAZ SCORES (Developmental Age Z-score)
# ============================================================================
# DAZ is a developmental assessment score from the GSED Long Form
# Extracted separately to allow for different analysis timepoints
daz <- gates %>%
  select(cohort, id, sex, age = age_daz, daz) %>%
  drop_na(daz) %>%
  arrange(cohort, id, age) %>%
  distinct(id, age, .keep_all = TRUE)

# ============================================================================
# COHORT-SPECIFIC PROCESSING
# ============================================================================

# Define brain regions for residual calculation
# These are the brain volume variables that will be adjusted for head circumference
brain_regions <- c("caudate", "putamen", "thalamus", "ventricles", 
                   "anterior_callosum", "central_callosum", "posterior_callosum", "corpus_callosum",
                   "supratentorial_csf", "supratentorial_tissue", 
                   "cerebellum", "globus_pallidus", 
                   "sc_csf", "sc_gray_matter", "sc_white_matter", "sc_corpus_callosum",
                   "sc_caudate", "sc_lentiform", "sc_hippocampus", "sc_thalamus")

# Define proper labels for visualization (currently not used but kept for reference)
brain_labels <- c(
  "Caudate", "Putamen", "Thalamus", "Ventricles",
  "Anterior Callosum", "Central Callosum", "Posterior Callosum", "Corpus Callosum",
  "Supratentorial CSF", "Supratentorial Tissue", 
  "Cerebellum", "Globus Pallidus",
  "SC CSF", "SC Gray Matter", "SC White Matter", "SC Corpus Callosum", 
  "SC Caudate", "SC Lentiform", "SC Hippocampus", "SC Thalamus"
)

# Select cohort for processing
# Uncomment the cohort you want to process:
cohorti <- "Khula"
# cohorti <- "SPACE_Explore"
# cohorti <- "PRIMES"
# cohorti <- "MINE"
# cohorti <- "PRISMA-Kenya"
# cohorti <- "PRISMA-Zambia"
# cohorti <- "PRISMA-Pakistan"

# Filter data for selected cohort and convert age from days to years
df <- gates %>%
  filter(cohort == cohorti) %>%
  mutate(age = age / 365.25)  # Convert age from days to years

# ============================================================================
# BRAIN VOLUME RESIDUAL CALCULATION
# ============================================================================
# Regress out head circumference (hcmri) effect from brain volume variables
# This creates head-size-adjusted brain volumes by removing the effect of 
# overall head size, which is important for comparing brain volumes across
# individuals with different head sizes
# 
# Method: Fit cubic polynomial model for each brain region, then calculate
# standardized residuals. Residuals represent brain volumes adjusted for head size.

for (var in brain_regions) {
  # Get complete cases for this specific variable pair (hcmri and brain region)
  complete_data <- df[!is.na(df$hcmri) & !is.na(df[[var]]), ]
  new_col_name <- paste0(var, "_residual")
  
  # Require at least 4 observations to fit cubic model (3 degrees of freedom + 1)
  if(nrow(complete_data) >= 4) {
    # Fit cubic polynomial model: brain_volume ~ poly(hcmri, 3)
    # Cubic allows for non-linear relationships between head size and brain volume
    model <- lm(as.formula(paste(var, "~ poly(hcmri, 3)")), data = complete_data)
    
    # Create a vector of NAs for the full dataset
    residuals_full <- rep(NA, nrow(df))
    
    # Fill in residuals only for complete cases
    residuals_full[!is.na(df$hcmri) & !is.na(df[[var]])] <- residuals(model)
    
    # Create new column with "_residual" suffix and standardized residuals
    # Standardization (z-scores) makes residuals comparable across brain regions
    df[[new_col_name]] <- scale(residuals_full)[,1]
  } else {
    # If insufficient data, set to NA
    df[[new_col_name]] <- NA
  }
}

# Save the processed cohort dataset with brain volume residuals
# Output location: data/processed/{cohort_name}.csv
write.csv(df, paste0("data/processed/", cohorti, ".csv"), row.names = FALSE)

# ============================================================================
# CREATE OUTCOME VARIABLES FOR ANALYSIS
# ============================================================================
# Process BSID and DAZ scores to create binary outcome variables (Normal vs Risk)
# Risk is defined as scores in the bottom 17th percentile (approximately 1 SD below mean)
# Only the latest assessment for each child is used

# BSID Cognitive: Convert to binary outcome (Normal/Risk) based on 17th percentile
bsidc0 <- bsidc %>% 
  mutate(age = age / 365.25) %>%  # Convert age from days to years
  filter(cohort == cohorti & !is.na(sex)) %>%
  # Keep only the latest measurement (highest age) for each repeated ID
  group_by(id) %>%
  filter(age == max(age)) %>%
  ungroup() %>%
  # Create binary outcome: Risk if score <= 17th percentile, Normal otherwise
  mutate(bsid = factor(ifelse(bsid > quantile(bsid, 0.17), 'Normal', 'Risk')))

# BSID Language: Convert to binary outcome
bsidl0 <- bsidl %>% 
  mutate(age = age / 365.25) %>%
  filter(cohort == cohorti & !is.na(sex)) %>%
  group_by(id) %>%
  filter(age == max(age)) %>%
  ungroup() %>%
  mutate(bsid = factor(ifelse(bsid > quantile(bsid, 0.17), 'Normal', 'Risk')))

# BSID Motor: Convert to binary outcome
bsidm0 <- bsidm %>% 
  mutate(age = age / 365.25) %>%
  filter(cohort == cohorti & !is.na(sex)) %>%
  group_by(id) %>%
  filter(age == max(age)) %>%
  ungroup() %>%
  mutate(bsid = factor(ifelse(bsid > quantile(bsid, 0.17), 'Normal', 'Risk')))

# DAZ: Convert to binary outcome
# Note: Only includes assessments at or after the threshold age (default: 1 year)
daz_second0 <- daz %>% 
  mutate(age = age / 365.25) %>%
  filter(cohort == cohorti & !is.na(sex) & age >= threshold) %>%
  group_by(id) %>%
  filter(age == max(age)) %>%
  ungroup() %>%
  mutate(daz = factor(ifelse(daz > quantile(daz, 0.17), 'Normal', 'Risk')))

# ============================================================================
# EXTRACT EARLY PREDICTORS
# ============================================================================
# Define predictor variables for early prediction models
# These include demographics, anthropometrics, brain volumes, and EEG features
vars <- c("age", "age_eeg", "medu", "ga", "fs", "ma", "hcmri", "bw", "bl", "waz", "haz", "whz",
          brain_regions, paste0(brain_regions, "_residual"), 
          "LowAlpha", "HighAlpha", "Beta", "Gamma", "ThetaBeta", "HighAlphaDelta", "AlphaPeak")

# Extract early predictor data (collected before threshold age, default: 1 year)
# This creates a dataset of early-life predictors that can be used to predict
# later developmental outcomes (BSID/DAZ scores)
early_predictors <- df %>%
  filter(age < threshold) %>%  # Only include measurements before threshold age
  select(id, all_of(vars)) %>%
  # For children with multiple early measurements, keep the one with most complete data
  group_by(id) %>%
  mutate(
    missing_count = rowSums(is.na(across(all_of(vars))))  # Count missing values in predictor columns
  ) %>%
  filter(missing_count == min(missing_count)) %>%  # Keep row with minimum missing values
  # If multiple rows have same missing count, keep the latest one (most recent measurement)
  filter(age == max(age)) %>%
  select(-missing_count) %>%
  ungroup() %>%
  mutate(age_eeg = age_eeg / 365.25)  # Convert EEG age from days to years

# ============================================================================
# MERGE OUTCOMES WITH EARLY PREDICTORS
# ============================================================================
# Combine outcome variables (BSID/DAZ scores) with early predictor variables
# This creates analysis-ready datasets where each row contains:
#   - Outcome: Latest BSID/DAZ assessment (binary: Normal/Risk)
#   - Predictors: Early-life measurements (before threshold age)
#   - Age information: Age at scan/assessment, age at outcome measurement

# Merge BSID Cognitive outcome with early predictors
bsidc <- bsidc0 %>%
  select(id, sex, age, bsid) %>%
  left_join(early_predictors, by = "id", suffix = c("_cog", "")) %>%
  relocate(age, .before = age_cog) %>%  # age = age at scan/early measurement
  relocate(age_eeg, .before = age_cog) %>%
  rename(age_scan = age) %>%  # Rename to clarify this is age at scan
  arrange(id)

# Merge BSID Language outcome with early predictors
bsidl <- bsidl0 %>%
  select(id, sex, age, bsid) %>%
  left_join(early_predictors, by = "id", suffix = c("_cog", "")) %>%
  relocate(age, .before = age_cog) %>%
  relocate(age_eeg, .before = age_cog) %>%
  rename(age_scan = age) %>%
  arrange(id)

# Merge BSID Motor outcome with early predictors
bsidm <- bsidm0 %>%
  select(id, sex, age, bsid) %>%
  left_join(early_predictors, by = "id", suffix = c("_cog", "")) %>%
  relocate(age, .before = age_cog) %>%
  relocate(age_eeg, .before = age_cog) %>%
  rename(age_scan = age) %>%
  arrange(id)

# Merge DAZ outcome with early predictors
daz_second <- daz_second0 %>%
  select(id, sex, age, daz) %>%
  left_join(early_predictors, by = "id", suffix = c("_cog", "")) %>%
  relocate(age, .before = age_cog) %>%
  relocate(age_eeg, .before = age_cog) %>%
  rename(age_scan = age) %>%
  arrange(id)

# ============================================================================
# SAVE FINAL PROCESSED DATASETS
# ============================================================================
# Save analysis-ready datasets to data/processed/ directory
# Each file contains outcome variable merged with early predictors
# Output files:
#   - {cohort}_bsidc.csv: BSID Cognitive outcome dataset
#   - {cohort}_bsidl.csv: BSID Language outcome dataset
#   - {cohort}_bsidm.csv: BSID Motor outcome dataset
#   - {cohort}_daz_second.csv: DAZ outcome dataset

write.csv(bsidc, paste0("data/processed/", cohorti, "_bsidc.csv"), row.names = FALSE)
write.csv(bsidl, paste0("data/processed/", cohorti, "_bsidl.csv"), row.names = FALSE)
write.csv(bsidm, paste0("data/processed/", cohorti, "_bsidm.csv"), row.names = FALSE)
write.csv(daz_second, paste0("data/processed/", cohorti, "_daz_second.csv"), row.names = FALSE)
