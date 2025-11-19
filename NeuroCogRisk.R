# ============================================================================
# Neurodevelopmental Risk Prediction Using XGBoost
# ============================================================================
# This script fits XGBoost models to predict neurodevelopmental outcomes
# (BSID scores and DAZ/GSED scores) using early-life predictors including:
#   - Demographics and birth characteristics
#   - Brain volume measurements (with head circumference adjustment)
#   - EEG frequency band features
#
# Models are fit separately for different cohorts and predictor combinations:
#   - PRIMES: BSID Cognitive, Language, Motor, and GSED
#   - Khula + SPACE_Explore: Combined cohort analysis
#   - MINE, PRISMA-Kenya, PRISMA-Zambia: GSED only
#   - All cohorts combined: GSED analysis
# ============================================================================

# Load required libraries
library(ggplot2)      # Plotting
library(patchwork)    # Combining plots
library(dplyr)        # Data manipulation
library(gridExtra)    # Table layout for plots
library(grid)         # Grid graphics for tables

# Source modeling functions from xgboost.R
# Note: This file contains fit_xgboost_miss(), create_shap_page_xgb(), and perf_row()
source("models/xgboost.R")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Function to save plots with 2x2 layout and performance table
# Creates a PDF with SHAP plots and performance summary table
save_plot_with_layout <- function(plot_list, performance_data, filename, 
                                  plot_title = "XGBoost SHAP Analysis", 
                                  table_title = "XGBoost Performance Summary", 
                                  width = 16, height = 10) {
  # Save to img/ directory (create if it doesn't exist)
  pdf(paste0("img/", filename), width = width, height = height)
  
  # Ensure we always have a 2x2 layout
  if (length(plot_list) == 1) {
    # Create empty plots for the remaining 3 positions
    empty_plot <- ggplot() + theme_void()
    plot_list <- c(plot_list, list(empty_plot, empty_plot, empty_plot))
  }
  
  print(wrap_plots(plot_list, ncol = 2, nrow = 2) + 
    plot_layout(guides = "collect") + theme(legend.position = "right") +
    plot_annotation(title = plot_title))
  
  tbl <- tableGrob(
    performance_data,
    rows = NULL, theme = ttheme_minimal(
      core = list(fg_params = list(cex = 0.9)),
      colhead = list(fg_params = list(fontface = "bold"))
    )
  )
  grid.arrange(
    tbl,
    top = textGrob(
      table_title,
      gp = gpar(fontsize = 16, fontface = "bold")
    )
  )
  dev.off()
}

# ============================================================================
# PREDICTOR VARIABLE DEFINITIONS
# ============================================================================
# Note on available outcomes by cohort:
#   - DAZ/GSED: Khula, SPACE_Explore, PRIMES, MINE, PRISMA-Kenya, PRISMA-Zambia
#   - BSID: Khula, PRIMES only

# Base demographic and anthropometric predictors
# These are available across most cohorts
vars_base <- c("sex", "medu", "ma", "ga", "bw", "bl", "waz", "haz", "whz", "age_scan")

# Brain regions from main field MRI (head circumference-adjusted residuals)
# MM = Main field measurements
regions <- c("caudate", "putamen", "thalamus", "ventricles", 
             "corpus_callosum",
             "supratentorial_csf", "supratentorial_tissue", 
             "cerebellum")

# Brain regions from spinal cord field MRI (head circumference-adjusted residuals)
# SC = Spinal cord field measurements
regions_sc <- c("sc_csf", "sc_gray_matter", "sc_white_matter", "sc_corpus_callosum",
                "sc_caudate", "sc_lentiform", "sc_hippocampus", "sc_thalamus")

# EEG frequency band features
eeg <- c("LowAlpha", "HighAlpha", "Beta", "Gamma", "ThetaBeta", "HighAlphaDelta", "AlphaPeak")

# ============================================================================
# PRIMES COHORT ANALYSIS
# ============================================================================
# PRIMES cohort has BSID Cognitive, Language, Motor, and GSED outcomes
# Note: Only 4 observations available for spinal cord (SC) volumes, so SC
#       regions are excluded from PRIMES analyses

cohorti <- "PRIMES"

# Load processed datasets from data/processed/ directory
bsidc <- read.csv(paste0("data/processed/", cohorti, "_bsidc.csv"))
bsidl <- read.csv(paste0("data/processed/", cohorti, "_bsidl.csv"))
bsidm <- read.csv(paste0("data/processed/", cohorti, "_bsidm.csv"))
daz <- read.csv(paste0("data/processed/", cohorti, "_daz_second.csv"))
# Prepare data: convert to factors and create binary maternal education variable
# Maternal education: Low (≤11 years) vs High (>11 years)
bsidc <- bsidc %>%
  mutate(sex = factor(sex), 
         medu = factor(ifelse(medu <= 11, 'Low', 'High'), levels = c('Low', 'High')),
         bsid = factor(bsid), 
         AlphaPeak = factor(AlphaPeak, levels = c(0, 1)))
bsidl <- bsidl %>%
  mutate(sex = factor(sex), 
         medu = factor(ifelse(medu <= 11, 'Low', 'High'), levels = c('Low', 'High')),
         bsid = factor(bsid), 
         AlphaPeak = factor(AlphaPeak, levels = c(0, 1)))
bsidm <- bsidm %>%
  mutate(sex = factor(sex), 
         medu = factor(ifelse(medu <= 11, 'Low', 'High'), levels = c('Low', 'High')),
         bsid = factor(bsid), 
         AlphaPeak = factor(AlphaPeak, levels = c(0, 1)))
daz <- daz %>%
  mutate(sex = factor(sex), 
         medu = factor(ifelse(medu <= 11, 'Low', 'High'), levels = c('Low', 'High')),
         daz = factor(daz), 
         AlphaPeak = factor(AlphaPeak, levels = c(0, 1)))

# Check sample sizes
c(nrow(bsidc), nrow(bsidl), nrow(bsidm), nrow(daz))

# Model 1: Baseline predictors only (excluding variables with high missingness)
# Remove ga (gestational age), haz (height-for-age Z-score), whz (weight-for-height Z-score)
vars1 <- setdiff(vars_base, c('ga', 'haz', 'whz'))

# Complete case analysis for baseline predictors
bsidc1 <- bsidc %>% drop_na(all_of(vars1))
bsidl1 <- bsidl %>% drop_na(all_of(vars1))
bsidm1 <- bsidm %>% drop_na(all_of(vars1))
daz1 <- daz %>% drop_na(all_of(vars1))
c(nrow(bsidc1), nrow(bsidl1), nrow(bsidm1), nrow(daz1))

# Fit XGBoost models with baseline predictors
# Note: Different seeds used for each outcome to ensure reproducibility
xgb_primes <- list()
xgb_primes$bsidc <- fit_xgboost_miss(bsidc1, "bsid", vars1, seed = 2)
xgb_primes$bsidl <- fit_xgboost_miss(bsidl1, "bsid", vars1, seed = 8)
xgb_primes$bsidm <- fit_xgboost_miss(bsidm1, "bsid", vars1, seed = 1)
xgb_primes$daz <- fit_xgboost_miss(daz1, "daz", setdiff(vars1, 'sex'), seed = 1)

# Create SHAP plots and performance summary
plot_xgb_primes <- create_shap_page_xgb(list(xgb_primes$bsidc, xgb_primes$bsidl, xgb_primes$bsidm, xgb_primes$daz), 
                                       c('BSID CCS', 'BSID LCS', 'BSID MCS', 'GSED'))
perf_xgb_primes <- dplyr::bind_rows(
  perf_row("BSID CCS", xgb_primes$bsidc),
  perf_row("BSID LCS", xgb_primes$bsidl),
  perf_row("BSID MCS", xgb_primes$bsidm),
  perf_row("GSED",      xgb_primes$daz)
)

# Save plots and performance table
save_plot_with_layout(
  plot_list = plot_xgb_primes,
  performance_data = perf_xgb_primes,
  filename = paste0(cohorti, ".pdf"),
  plot_title = "XGBoost SHAP Analysis - PRIMES",
  table_title = "XGBoost Performance Summary - PRIMES"
)

# Model 2: Baseline + Brain volume residuals (MM = Main field MRI)
# Use head circumference-adjusted brain volume residuals
vars_residual <- c(vars_base, paste0(regions, "_residual"))
vars_residual <- setdiff(vars_residual, c('ga', 'haz', 'whz'))

# Complete case analysis for brain volume residuals
bsidc1 <- bsidc %>% drop_na(all_of(paste0(regions, "_residual")))
bsidl1 <- bsidl %>% drop_na(all_of(paste0(regions, "_residual")))
bsidm1 <- bsidm %>% drop_na(all_of(paste0(regions, "_residual")))
daz1 <- daz %>% drop_na(all_of(paste0(regions, "_residual")))
c(nrow(bsidc1), nrow(bsidl1), nrow(bsidm1), nrow(daz1))

# Fit XGBoost models with baseline + brain volume predictors
xgb_primes_mm <- list()
xgb_primes_mm$bsidc <- fit_xgboost_miss(bsidc1, "bsid", vars_residual, seed = 2)
xgb_primes_mm$bsidl <- fit_xgboost_miss(bsidl1, "bsid", vars_residual, seed = 8)
xgb_primes_mm$bsidm <- fit_xgboost_miss(bsidm1, "bsid", vars_residual, seed = 1)
xgb_primes_mm$daz <- fit_xgboost_miss(daz1, "daz", setdiff(vars_residual, 'sex'), seed = 1)

plot_xgb_primes_mm <- create_shap_page_xgb(list(xgb_primes_mm$bsidc, xgb_primes_mm$bsidl, xgb_primes_mm$bsidm, xgb_primes_mm$daz), c('BSID CCS', 'BSID LCS', 'BSID MCS', 'GSED'))
perf_xgb_primes_mm <- dplyr::bind_rows(
  perf_row("BSID CCS", xgb_primes_mm$bsidc),
  perf_row("BSID LCS", xgb_primes_mm$bsidl),
  perf_row("BSID MCS", xgb_primes_mm$bsidm),
  perf_row("GSED",      xgb_primes_mm$daz)
)

save_plot_with_layout(
  plot_list = plot_xgb_primes_mm,
  performance_data = perf_xgb_primes_mm,
  filename = paste0(cohorti, "_MM.pdf"),
  plot_title = "XGBoost SHAP Analysis - PRIMES - MM",
  table_title = "XGBoost Performance Summary - PRIMES - MM"
)

# ============================================================================
# KHULA + SPACE_Explore COMBINED COHORT ANALYSIS
# ============================================================================
# Combined analysis of Khula and SPACE_Explore cohorts
# Khula has BSID outcomes; both cohorts have GSED outcomes

cohorti <- "Khula+SPACE"

# Load Khula BSID datasets
bsidc <- read.csv("data/processed/Khula_bsidc.csv")
bsidl <- read.csv("data/processed/Khula_bsidl.csv")
bsidm <- read.csv("data/processed/Khula_bsidm.csv")

# Combine GSED/DAZ from both cohorts
daz <- rbind(read.csv("data/processed/Khula_daz_second.csv"), 
             read.csv("data/processed/SPACE_Explore_daz_second.csv"))
bsidc <- bsidc %>%
  mutate(sex = factor(sex),
         medu = factor(ifelse(medu <= 11, 'Low', 'High'), levels = c('Low', 'High')),
         bsid = factor(bsid), 
         AlphaPeak = factor(AlphaPeak, levels = c(0, 1)))
bsidl <- bsidl %>%
  mutate(sex = factor(sex),
         medu = factor(ifelse(medu <= 11, 'Low', 'High'), levels = c('Low', 'High')), 
         bsid = factor(bsid), 
         AlphaPeak = factor(AlphaPeak, levels = c(0, 1)))
bsidm <- bsidm %>%
  mutate(sex = factor(sex),
         medu = factor(ifelse(medu <= 11, 'Low', 'High'), levels = c('Low', 'High')), 
         bsid = factor(bsid), 
         AlphaPeak = factor(AlphaPeak, levels = c(0, 1)))
daz <- daz %>%
  mutate(sex = factor(sex), 
         medu = factor(ifelse(medu <= 11, 'Low', 'High'), levels = c('Low', 'High')), 
         daz = factor(daz), 
         AlphaPeak = factor(AlphaPeak, levels = c(0, 1)))
c(nrow(bsidc), nrow(bsidl), nrow(bsidm), nrow(daz))
summary(daz$age_cog)

# Model 1: Baseline predictors only
# Complete case analysis for all baseline variables
bsidc1 <- bsidc %>% drop_na(all_of(vars_base))
bsidl1 <- bsidl %>% drop_na(all_of(vars_base))
bsidm1 <- bsidm %>% drop_na(all_of(vars_base))
daz1 <- daz %>% drop_na(all_of(vars_base))
c(nrow(bsidc1), nrow(bsidl1), nrow(bsidm1), nrow(daz1))

# Fit XGBoost models with baseline predictors
# Note: Sex excluded from GSED model (may have been problematic in this cohort)
xgb_uct <- list()
xgb_uct$bsidc <- fit_xgboost_miss(bsidc1, "bsid", vars_base)
xgb_uct$bsidl <- fit_xgboost_miss(bsidl1, "bsid", vars_base)
xgb_uct$bsidm <- fit_xgboost_miss(bsidm1, "bsid", vars_base)
xgb_uct$daz <- fit_xgboost_miss(daz1, "daz", setdiff(vars_base, 'sex'))

plot_xgb_uct <- create_shap_page_xgb(list(xgb_uct$bsidc, xgb_uct$bsidl, xgb_uct$bsidm, xgb_uct$daz), c('BSID CCS', 'BSID LCS', 'BSID MCS', 'GSED'))
perf_xgb_uct <- dplyr::bind_rows(
  perf_row("BSID CCS", xgb_uct$bsidc),
  perf_row("BSID LCS", xgb_uct$bsidl),
  perf_row("BSID MCS", xgb_uct$bsidm),
  perf_row("GSED",      xgb_uct$daz)
)

save_plot_with_layout(
  plot_list = plot_xgb_uct,
  performance_data = perf_xgb_uct,
  filename = paste0(cohorti, ".pdf"),
  plot_title = "XGBoost SHAP Analysis - Khula + SPACE",
  table_title = "XGBoost Performance Summary - Khula + SPACE"
)

# Model 2: Baseline + EEG features
# Combine baseline predictors with EEG frequency band features
vars_residual <- c(vars_base, eeg)

# Complete case analysis for EEG variables
bsidc1 <- bsidc %>% drop_na(all_of(eeg))
bsidl1 <- bsidl %>% drop_na(all_of(eeg))
bsidm1 <- bsidm %>% drop_na(all_of(eeg))
daz1 <- daz %>% drop_na(all_of(eeg))
c(nrow(bsidc1), nrow(bsidl1), nrow(bsidm1), nrow(daz1))

# Fit XGBoost models with baseline + EEG predictors
# Increased max_trials for more thorough parameter search
xgb_uct_eeg <- list()
xgb_uct_eeg$bsidc <- fit_xgboost_miss(bsidc1, "bsid", vars_residual, max_trials = 2000)
xgb_uct_eeg$bsidl <- fit_xgboost_miss(bsidl1, "bsid", vars_residual, max_trials = 2000)
xgb_uct_eeg$bsidm <- fit_xgboost_miss(bsidm1, "bsid", vars_residual, max_trials = 2000)
xgb_uct_eeg$daz <- fit_xgboost_miss(daz1, "daz", setdiff(vars_residual, 'sex'), max_trials = 2000)

plot_xgb_uct_eeg <- create_shap_page_xgb(list(xgb_uct_eeg$bsidc, xgb_uct_eeg$bsidl, xgb_uct_eeg$bsidm, xgb_uct_eeg$daz), c('BSID CCS', 'BSID LCS', 'BSID MCS', 'GSED'))
perf_xgb_uct_eeg <- dplyr::bind_rows(
  perf_row("BSID CCS", xgb_uct_eeg$bsidc),
  perf_row("BSID LCS", xgb_uct_eeg$bsidl),
  perf_row("BSID MCS", xgb_uct_eeg$bsidm),
  perf_row("GSED",      xgb_uct_eeg$daz)
)

save_plot_with_layout(
  plot_list = plot_xgb_uct_eeg,
  performance_data = perf_xgb_uct_eeg,
  filename = paste0(cohorti, "_EEG.pdf"),
  plot_title = "XGBoost SHAP Analysis - Khula + SPACE - EEG",
  table_title = "XGBoost Performance Summary - Khula + SPACE - EEG"
)

# Model 3: Baseline + Brain volume residuals (MM + SC)
# Combine baseline predictors with both main field and spinal cord brain volume residuals
vars_residual <- c(vars_base, paste0(c(regions, regions_sc), "_residual"))

# Complete case analysis for brain volume residuals
bsidc1 <- bsidc %>% drop_na(all_of(paste0(c(regions, regions_sc), "_residual")))
bsidl1 <- bsidl %>% drop_na(all_of(paste0(c(regions, regions_sc), "_residual")))
bsidm1 <- bsidm %>% drop_na(all_of(paste0(c(regions, regions_sc), "_residual")))
daz1 <- daz %>% drop_na(all_of(paste0(c(regions, regions_sc), "_residual")))
c(nrow(bsidc1), nrow(bsidl1), nrow(bsidm1), nrow(daz1))

# Fit XGBoost models with baseline + MRI predictors
xgb_uct_mri <- list()
xgb_uct_mri$bsidc <- fit_xgboost_miss(bsidc1, "bsid", vars_residual, seed = 1)
xgb_uct_mri$bsidl <- fit_xgboost_miss(bsidl1, "bsid", vars_residual, max_trials = 2000, seed = 3)
xgb_uct_mri$bsidm <- fit_xgboost_miss(bsidm1, "bsid", vars_residual, max_trials = 2000)
xgb_uct_mri$daz <- fit_xgboost_miss(daz1, "daz", setdiff(vars_residual, 'sex'), max_trials = 2000, seed = 2)

plot_xgb_uct_mri <- create_shap_page_xgb(list(xgb_uct_mri$bsidc, xgb_uct_mri$bsidl, xgb_uct_mri$bsidm, xgb_uct_mri$daz), c('BSID CCS', 'BSID LCS', 'BSID MCS', 'GSED'))
perf_xgb_uct_mri <- dplyr::bind_rows(
  perf_row("BSID CCS", xgb_uct_mri$bsidc),
  perf_row("BSID LCS", xgb_uct_mri$bsidl),
  perf_row("BSID MCS", xgb_uct_mri$bsidm),
  perf_row("GSED",      xgb_uct_mri$daz)
)

save_plot_with_layout(
  plot_list = plot_xgb_uct_mri,
  performance_data = perf_xgb_uct_mri,
  filename = paste0(cohorti, "_MM_SC.pdf"),
  plot_title = "XGBoost SHAP Analysis - Khula + SPACE - MRI",
  table_title = "XGBoost Performance Summary - Khula + SPACE - MRI"
)

# Model 4: Baseline + EEG + Brain volume residuals (MM + SC)
# Full model with all available predictors
vars_residual <- c(vars_base, eeg, paste0(c(regions, regions_sc), "_residual"))

# Complete case analysis for all predictors (EEG + MRI)
bsidc1 <- bsidc %>% drop_na(all_of(c(eeg, paste0(c(regions, regions_sc), "_residual"))))
bsidl1 <- bsidl %>% drop_na(all_of(c(eeg, paste0(c(regions, regions_sc), "_residual"))))
bsidm1 <- bsidm %>% drop_na(all_of(c(eeg, paste0(c(regions, regions_sc), "_residual"))))
daz1 <- daz %>% drop_na(all_of(c(eeg, paste0(c(regions, regions_sc), "_residual"))))
c(nrow(bsidc1), nrow(bsidl1), nrow(bsidm1), nrow(daz1))

# Fit XGBoost models with all predictors
xgb_uct_mri_eeg <- list()
xgb_uct_mri_eeg$bsidc <- fit_xgboost_miss(bsidc1, "bsid", vars_residual, max_trials = 2000, seed = 1)
xgb_uct_mri_eeg$bsidl <- fit_xgboost_miss(bsidl1, "bsid", vars_residual, max_trials = 2000, seed = 2)
xgb_uct_mri_eeg$bsidm <- fit_xgboost_miss(bsidm1, "bsid", vars_residual, max_trials = 2000)
xgb_uct_mri_eeg$daz <- fit_xgboost_miss(daz1, "daz", setdiff(vars_residual, 'sex'), max_trials = 2000, seed = 1)

plot_xgb_uct_mri_eeg <- create_shap_page_xgb(list(xgb_uct_mri_eeg$bsidc, xgb_uct_mri_eeg$bsidl, xgb_uct_mri_eeg$bsidm, xgb_uct_mri_eeg$daz), c('BSID CCS', 'BSID LCS', 'BSID MCS', 'GSED'))
perf_xgb_uct_mri_eeg <- dplyr::bind_rows(
  perf_row("BSID CCS", xgb_uct_mri_eeg$bsidc),
  perf_row("BSID LCS", xgb_uct_mri_eeg$bsidl),
  perf_row("BSID MCS", xgb_uct_mri_eeg$bsidm),
  perf_row("GSED",      xgb_uct_mri_eeg$daz)
)

save_plot_with_layout(
  plot_list = plot_xgb_uct_mri_eeg,
  performance_data = perf_xgb_uct_mri_eeg,
  filename = paste0(cohorti, "_EEG_MM_SC.pdf"),
  plot_title = "XGBoost SHAP Analysis - Khula + SPACE - EEG + MRI",
  table_title = "XGBoost Performance Summary - Khula + SPACE - EEG + MRI"
)

# Model 5: Baseline + Brain volume residuals (MM only, no SC)
# Main field MRI brain volumes only (excluding spinal cord measurements)
vars_residual <- c(vars_base, paste0(regions, "_residual"))

# Complete case analysis for main field brain volume residuals
bsidc1 <- bsidc %>% drop_na(all_of(paste0(regions, "_residual")))
bsidl1 <- bsidl %>% drop_na(all_of(paste0(regions, "_residual")))
bsidm1 <- bsidm %>% drop_na(all_of(paste0(regions, "_residual")))
daz1 <- daz %>% drop_na(all_of(paste0(regions, "_residual")))

# Fit XGBoost models with baseline + main field MRI predictors
xgb_uct_mm <- list()
xgb_uct_mm$bsidc <- fit_xgboost_miss(bsidc1, "bsid", vars_residual, max_trials = 2000)
xgb_uct_mm$bsidl <- fit_xgboost_miss(bsidl1, "bsid", vars_residual, max_trials = 2000)
xgb_uct_mm$bsidm <- fit_xgboost_miss(bsidm1, "bsid", vars_residual, max_trials = 2000)
xgb_uct_mm$daz <- fit_xgboost_miss(daz1, "daz", vars_residual, max_trials = 2000)

plot_xgb_uct_mm <- create_shap_page_xgb(list(xgb_uct_mm$bsidc, xgb_uct_mm$bsidl, xgb_uct_mm$bsidm, xgb_uct_mm$daz), c('BSID CCS', 'BSID LCS', 'BSID MCS', 'GSED'))
perf_xgb_uct_mm <- dplyr::bind_rows(
  perf_row("BSID CCS", xgb_uct_mm$bsidc),
  perf_row("BSID LCS", xgb_uct_mm$bsidl),
  perf_row("BSID MCS", xgb_uct_mm$bsidm),
  perf_row("GSED",      xgb_uct_mm$daz)
)

save_plot_with_layout(
  plot_list = plot_xgb_uct_mm,
  performance_data = perf_xgb_uct_mm,
  filename = paste0(cohorti, "_MM.pdf"),
  plot_title = "XGBoost SHAP Analysis - Khula + SPACE - MM",
  table_title = "XGBoost Performance Summary - Khula + SPACE - MM"
)

# Model 6: Baseline + EEG + Brain volume residuals (MM only)
# Main field MRI brain volumes combined with EEG features
vars_residual <- c(vars_base, eeg, paste0(regions, "_residual"))

# Complete case analysis for EEG + main field brain volume residuals
bsidc1 <- bsidc %>% drop_na(all_of(c(eeg, paste0(regions, "_residual"))))
bsidl1 <- bsidl %>% drop_na(all_of(c(eeg, paste0(regions, "_residual"))))
bsidm1 <- bsidm %>% drop_na(all_of(c(eeg, paste0(regions, "_residual"))))
daz1 <- daz %>% drop_na(all_of(c(eeg, paste0(regions, "_residual"))))
c(nrow(bsidc1), nrow(bsidl1), nrow(bsidm1), nrow(daz1))

# Fit XGBoost models with baseline + EEG + main field MRI predictors
xgb_uct_mm_eeg <- list()
xgb_uct_mm_eeg$bsidc <- fit_xgboost_miss(bsidc1, "bsid", vars_residual, max_trials = 2000, seed = 2)
xgb_uct_mm_eeg$bsidl <- fit_xgboost_miss(bsidl1, "bsid", vars_residual, max_trials = 2000, seed = 1)
xgb_uct_mm_eeg$bsidm <- fit_xgboost_miss(bsidm1, "bsid", vars_residual, max_trials = 2000)
xgb_uct_mm_eeg$daz <- fit_xgboost_miss(daz1, "daz", setdiff(vars_residual, 'sex'), max_trials = 2000)

plot_xgb_uct_mm_eeg <- create_shap_page_xgb(list(xgb_uct_mm_eeg$bsidc, xgb_uct_mm_eeg$bsidl, xgb_uct_mm_eeg$bsidm, xgb_uct_mm_eeg$daz), c('BSID CCS', 'BSID LCS', 'BSID MCS', 'GSED'))
perf_xgb_uct_mm_eeg <- dplyr::bind_rows(
  perf_row("BSID CCS", xgb_uct_mm_eeg$bsidc),
  perf_row("BSID LCS", xgb_uct_mm_eeg$bsidl),
  perf_row("BSID MCS", xgb_uct_mm_eeg$bsidm),
  perf_row("GSED",      xgb_uct_mm_eeg$daz)
)

save_plot_with_layout(
  plot_list = plot_xgb_uct_mm_eeg,
  performance_data = perf_xgb_uct_mm_eeg,
  filename = paste0(cohorti, "_EEG_MM.pdf"),
  plot_title = "XGBoost SHAP Analysis - Khula + SPACE - EEG + MM",
  table_title = "XGBoost Performance Summary - Khula + SPACE - EEG + MM"
)

# ============================================================================
# MINE COHORT ANALYSIS
# ============================================================================
# MINE cohort has GSED/DAZ outcome only

cohorti <- "MINE"
daz <- read.csv(paste0("data/processed/", cohorti, "_daz_second.csv"))
daz <- daz %>%
  mutate(sex = factor(sex), 
         medu = factor(ifelse(medu <= 11, 'Low', 'High'), levels = c('Low', 'High')),
         daz = factor(daz), 
         AlphaPeak = factor(AlphaPeak, levels = c(0, 1)))
nrow(daz)

daz1 <- daz %>% drop_na(all_of(vars_base))
nrow(daz1)

xgb_mine <- list()
xgb_mine$daz <- fit_xgboost_miss(daz1, "daz", setdiff(vars_base, 'sex'))

plot_xgb_mine <- create_shap_page_xgb(xgb_mine$daz, 'GSED')
perf_xgb_mine <- perf_row("GSED", xgb_mine$daz)

save_plot_with_layout(
  plot_list = plot_xgb_mine,
  performance_data = perf_xgb_mine,
  filename = paste0(cohorti, ".pdf"),
  plot_title = "XGBoost SHAP Analysis - MINE",
  table_title = "XGBoost Performance Summary - MINE"
)

vars_residual <- c(vars_base, paste0(c(regions, regions_sc), "_residual"))
daz1 <- daz %>% drop_na(all_of(paste0(c(regions, regions_sc), "_residual")))
nrow(daz1)

xgb_mine_mri <- list()
xgb_mine_mri$daz <- fit_xgboost_miss(daz1, "daz", setdiff(vars_residual, 'sex'), seed = 1)

plot_xgb_mine_mri <- create_shap_page_xgb(xgb_mine_mri$daz, "GSED")
perf_xgb_mine_mri <- perf_row("GSED", xgb_mine_mri$daz)

save_plot_with_layout(
  plot_list = plot_xgb_mine_mri,
  performance_data = perf_xgb_mine_mri,
  filename = paste0(cohorti, "_MM_SC.pdf"),
  plot_title = "XGBoost SHAP Analysis - MINE - MRI",
  table_title = "XGBoost Performance Summary - MINE - MRI"
)

vars_residual <- c(vars_base, paste0(regions, "_residual"))
daz1 <- daz %>% drop_na(all_of(paste0(regions, "_residual")))
nrow(daz1)

xgb_mine_mm <- list()
xgb_mine_mm$daz <- fit_xgboost_miss(daz1, "daz", setdiff(vars_residual, 'sex'), max_trials = 2000, seed = 1)

plot_xgb_mine_mm <- create_shap_page_xgb(xgb_mine_mm$daz, "GSED")
perf_xgb_mine_mm <- perf_row("GSED", xgb_mine_mm$daz)

save_plot_with_layout(
  plot_list = plot_xgb_mine_mm,
  performance_data = perf_xgb_mine_mm,
  filename = paste0(cohorti, "_MM.pdf"),
  plot_title = "XGBoost SHAP Analysis - MINE - MM",
  table_title = "XGBoost Performance Summary - MINE - MM"
)

vars_residual <- c(vars_base, paste0(regions_sc, "_residual"))
daz1 <- daz %>% drop_na(all_of(paste0(regions_sc, "_residual")))
nrow(daz1)

xgb_mine_sc <- list()
xgb_mine_sc$daz <- fit_xgboost_miss(daz1, "daz", setdiff(vars_residual, 'sex'), max_trials = 2000, seed = 1)

plot_xgb_mine_sc <- create_shap_page_xgb(xgb_mine_sc$daz, "GSED")
perf_xgb_mine_sc <- perf_row("GSED", xgb_mine_sc$daz)

save_plot_with_layout(
  plot_list = plot_xgb_mine_sc,
  performance_data = perf_xgb_mine_sc,
  filename = paste0(cohorti, "_SC.pdf"),
  plot_title = "XGBoost SHAP Analysis - MINE - SC",
  table_title = "XGBoost Performance Summary - MINE - SC"
)

# ============================================================================
# PRISMA-ZAMBIA COHORT ANALYSIS
# ============================================================================
# PRISMA-Zambia cohort has GSED/DAZ outcome only

cohorti <- "PRISMA-Zambia"
daz <- read.csv(paste0("data/processed/", cohorti, "_daz_second.csv"))
daz <- daz %>%
  mutate(sex = factor(sex), 
         medu = factor(ifelse(medu <= 11, 'Low', 'High'), levels = c('Low', 'High')),
         daz = factor(daz), 
         AlphaPeak = factor(AlphaPeak, levels = c(0, 1)))
nrow(daz)

vars1 <- setdiff(vars_base, c('waz', 'haz', 'whz'))
daz1 <- daz %>% drop_na(all_of(vars1))
nrow(daz1)

xgb_zambia <- list()
xgb_zambia$daz <- fit_xgboost_miss(daz1, "daz", setdiff(vars1, 'sex'))

plot_xgb_zambia <- create_shap_page_xgb(xgb_zambia$daz, 'GSED')
perf_xgb_zambia <- perf_row("GSED", xgb_zambia$daz)

save_plot_with_layout(
  plot_list = plot_xgb_zambia,
  performance_data = perf_xgb_zambia,
  filename = paste0(cohorti, ".pdf"),
  plot_title = "XGBoost SHAP Analysis - Zambia",
  table_title = "XGBoost Performance Summary - Zambia"
)

vars_residual <- c(vars_base, paste0(c(regions, regions_sc), "_residual"))
daz1 <- daz %>% drop_na(all_of(paste0(c(regions, regions_sc), "_residual")))
nrow(daz1)

xgb_zambia_mri <- list()
xgb_zambia_mri$daz <- fit_xgboost_miss(daz1, "daz", 
                                       setdiff(vars_residual, c('sex', 'waz', 'haz', 'whz')), 
                                       max_trials = 2000)

plot_xgb_zambia_mri <- create_shap_page_xgb(xgb_zambia_mri$daz, "GSED")
perf_xgb_zambia_mri <- perf_row("GSED", xgb_zambia_mri$daz)

save_plot_with_layout(
  plot_list = plot_xgb_zambia_mri,
  performance_data = perf_xgb_zambia_mri,
  filename = paste0(cohorti, "_MM_SC.pdf"),
  plot_title = "XGBoost SHAP Analysis - PRISMA Zambia - MRI",
  table_title = "XGBoost Performance Summary - PRISMA Zambia - MRI"
)

vars_residual <- c(vars_base, paste0(regions, "_residual"))
daz1 <- daz %>% drop_na(all_of(paste0(regions, "_residual")))
nrow(daz1)

xgb_zambia_mm <- list()
xgb_zambia_mm$daz <- fit_xgboost_miss(daz1, "daz", 
                                      setdiff(vars_residual, c('sex', 'waz', 'haz', 'whz')), 
                                      max_trials = 2000)

plot_xgb_zambia_mm <- create_shap_page_xgb(xgb_zambia_mm$daz, "GSED")
perf_xgb_zambia_mm <- perf_row("GSED", xgb_zambia_mm$daz)

save_plot_with_layout(
  plot_list = plot_xgb_zambia_mm,
  performance_data = perf_xgb_zambia_mm,
  filename = paste0(cohorti, "_MM.pdf"),
  plot_title = "XGBoost SHAP Analysis - PRISMA Zambia - MM",
  table_title = "XGBoost Performance Summary - PRISMA Zambia - MM"
)

vars_residual <- c(vars_base, paste0(regions_sc, "_residual"))
daz1 <- daz %>% drop_na(all_of(paste0(regions_sc, "_residual")))
nrow(daz1)

xgb_zambia_sc <- list()
xgb_zambia_sc$daz <- fit_xgboost_miss(daz1, "daz", 
                                      setdiff(vars_residual, c('sex', 'waz', 'haz', 'whz')), 
                                      max_trials = 2000)

plot_xgb_zambia_sc <- create_shap_page_xgb(xgb_zambia_sc$daz, "GSED")
perf_xgb_zambia_sc <- perf_row("GSED", xgb_zambia_sc$daz)

save_plot_with_layout(
  plot_list = plot_xgb_zambia_sc,
  performance_data = perf_xgb_zambia_sc,
  filename = paste0(cohorti, "_SC.pdf"),
  plot_title = "XGBoost SHAP Analysis - PRISMA Zambia - SC",
  table_title = "XGBoost Performance Summary - PRISMA Zambia - SC"
)

# ============================================================================
# PRISMA-KENYA COHORT ANALYSIS
# ============================================================================
# PRISMA-Kenya cohort has GSED/DAZ outcome only

cohorti <- "PRISMA-Kenya"
daz <- read.csv(paste0("data/processed/", cohorti, "_daz_second.csv"))
daz <- daz %>%
  mutate(sex = factor(sex), 
         medu = factor(ifelse(medu <= 11, 'Low', 'High'), levels = c('Low', 'High')),
         daz = factor(daz), 
         AlphaPeak = factor(AlphaPeak, levels = c(0, 1)))
nrow(daz)

vars1 <- setdiff(vars_base, c('waz', 'haz', 'whz'))
daz1 <- daz %>% drop_na(all_of(vars1))
nrow(daz1)

xgb_kenya <- list()
xgb_kenya$daz <- fit_xgboost_miss(daz1, "daz", setdiff(vars1, 'sex'))

plot_xgb_kenya <- create_shap_page_xgb(xgb_kenya$daz, 'GSED')
perf_xgb_kenya <- perf_row("GSED", xgb_kenya$daz)

save_plot_with_layout(
  plot_list = plot_xgb_kenya,
  performance_data = perf_xgb_kenya,
  filename = paste0(cohorti, ".pdf"),
  plot_title = "XGBoost SHAP Analysis - Kenya",
  table_title = "XGBoost Performance Summary - Kenya"
)

vars_residual <- c(vars_base, paste0(regions, "_residual"))
daz1 <- daz %>% drop_na(all_of(paste0(regions, "_residual")))
nrow(daz1)

xgb_kenya_mm <- list()
xgb_kenya_mm$daz <- fit_xgboost_miss(daz1, "daz", 
                                     setdiff(vars_residual, c('sex', 'waz', 'haz', 'whz')), 
                                     max_trials = 2000, seed = 1)

plot_xgb_kenya_mm <- create_shap_page_xgb(xgb_kenya_mm$daz, "GSED")
perf_xgb_kenya_mm <- perf_row("GSED", xgb_kenya_mm$daz)

save_plot_with_layout(
  plot_list = plot_xgb_kenya_mm,
  performance_data = perf_xgb_kenya_mm,
  filename = paste0(cohorti, "_MM.pdf"),
  plot_title = "XGBoost SHAP Analysis - PRISMA Kenya",
  table_title = "XGBoost Performance Summary - PRISMA Kenya"
)
