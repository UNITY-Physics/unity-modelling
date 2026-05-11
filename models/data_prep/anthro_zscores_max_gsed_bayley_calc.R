##########################################################
## Standardizing cohort data for analyses - for FlyWheel upload
## Datasets that are already in FlyWheel template structure
## Brain Health Metrics IHME
## May 2026
##########################################################

## Set up ---------------------------------------------------------------------------------------------------
rm(list=ls())
username <- Sys.getenv("USER")

library(readxl)
library(tidyr)
library(ggplot2)
library(readstata13)
library(data.table)
library(dplyr)
library(survey)
library(haven)
library(anthro, lib.loc = paste0("/FILEPATH/", username, "/"))
library(anthroplus, lib.loc = paste0("/FILEPATH/", username, "/"))

#System configuration
l_root <- "/FILEPATH/"
j_root <- "/FILEPATH/"

#Set objects
file_path <- "/FILEPATH/"
file_name <- "UCT-Khula-Hyperfine_MERGED_EEG_sc.csv"

save_dir <- paste0(l_root, "FILEPATH/")

code_dir <- paste0("/FILEPATH/", username, "/FILEPATH/")
source(paste0(code_dir, "process_cohort_data_unity/standardize_cohort_data_functions.R"))

## Read in data ----------------------------------------------------------------------
dt_orig <- as.data.table(read.csv(paste0(file_path, file_name)))
#dt_orig <- dt_orig[!CohortName == "", ]

dt_orig[, UniqueStudyID := sub("_[^_]+$", "", StudyID)]
dt_orig[childTimepointHC_MRI_cm==0, childTimepointHC_MRI_cm := NA]

## Calculate maternal education in continuous years from HollingsHead scale ----------------------------------------------------------------------
#RESONANCE
# dt_orig[, MaternalEducation_schoolingYears := fcase(
#   MaternalEducation_HHS == 1, 6,
#   MaternalEducation_HHS == 2, 9,
#   MaternalEducation_HHS == 3, 11,
#   MaternalEducation_HHS == 4, 12,
#   MaternalEducation_HHS == 5, 14,
#   MaternalEducation_HHS == 6, 16,
#   MaternalEducation_HHS == 7, 18,
#   default = NA_real_
# )]

#CEL
# dt_orig[, MaternalEducation_schoolingYears_old := MaternalEducation_schoolingYears]
# dt_orig[, MaternalEducation_schoolingYears := NA]
# dt_orig[, MaternalEducation_schoolingYears := as.numeric(MaternalEducation_schoolingYears)]
# dt_orig[, MaternalEducation_schoolingYears := fcase(
#   MaternalEducation_schoolingYears_old == "", NA_real_,
#   MaternalEducation_schoolingYears_old == "illiterate", 0,
#   MaternalEducation_schoolingYears_old == "below_primary", 3,
#   MaternalEducation_schoolingYears_old == "primary_pass", 5,
#   MaternalEducation_schoolingYears_old == "highschool_pass", 10,
#   MaternalEducation_schoolingYears_old == "inter_pass", 12,
#   MaternalEducation_schoolingYears_old == "polytechnic_or_diploma", 14,
#   MaternalEducation_schoolingYears_old == "graduate", 15,
#   MaternalEducation_schoolingYears_old == "postgraduate", 17,
#   default = NA_real_
# )]
# 
# dt_orig <- subset(dt_orig, select = -c(MaternalEducation_schoolingYears_old))

#BCD and ENAT
# dt_orig[, MaternalEducation_schoolingYears_old := MaternalEducation_schoolingYears]
# dt_orig[, MaternalEducation_schoolingYears := NA]
# dt_orig[, MaternalEducation_schoolingYears := as.numeric(MaternalEducation_schoolingYears)]
# dt_orig[, MaternalEducation_schoolingYears := fcase(
#   MaternalEducation_schoolingYears_old =="No education", 0,
#   MaternalEducation_schoolingYears_old == "Able to read and write", 3,
#   MaternalEducation_schoolingYears_old == "Primary", 4,
#   MaternalEducation_schoolingYears_old == "Secondary", 10,
#   MaternalEducation_schoolingYears_old == "Technical/vocational", 12,
#   MaternalEducation_schoolingYears_old == "Higher", 15,
#   default = NA_real_
# )]
# 
# dt_orig <- subset(dt_orig, select = -c(MaternalEducation_schoolingYears_old))

#PRISMA Kintampo
table(dt_orig$MaternalEducation_schoolingYears)
dt_orig[MaternalEducation_schoolingYears==77, MaternalEducation_schoolingYears := 0]
table(dt_orig$MaternalEducation_schoolingYears)

## Calculate anthro z-scores and binary variables ----------------------------------------------------------------------
process_anthro <- function(df) {
  
  # 0. Calculate raw BMI directly safely (ALWAYS add the column)
  if (all(c("weight_kgs", "height_cm") %in% names(df))) {
    df[, bmi := as.numeric(weight_kgs) / ((as.numeric(height_cm) / 100)^2)]
    
    if (all(is.na(df$bmi))) {
      message("Notice: Calculated 'bmi' is entirely NA, but column was still added.")
    }
  } else {
    message("Notice: 'weight_kgs' or 'height_cm' missing. 'bmi' column added as all NA.")
    df[, bmi := NA_real_]
  }
  
  # 1. Check for absolute minimum required variables
  if (!all(c("child_sex", "other_age_days") %in% names(df))) {
    message("Notice: Missing 'child_sex' or 'other_age_days'. Skipping z-score calculations entirely.")
    return(df)
  }
  
  # Only run if output columns don't already exist
  if (!all(c("zlen", "zwei", "zwfl", "zbmi", "zhc") %in% names(df))) {
    
    sex_vec <- fcase(
      df$child_sex == "Male", 1L,
      df$child_sex == "Female", 2L,
      default = NA_integer_
    )
    
    # Keep the true age for the dataset, but create a safe calculation age
    age_days <- as.numeric(df$other_age_days)
    calc_age_days <- copy(age_days)
    
    if (all(is.na(sex_vec)) || all(is.na(age_days))) {
      message("Notice: All rows have missing sex or age data. Skipping z-score calculations.")
    } else {
      
      weight_kg <- if ("weight_kgs" %in% names(df)) as.numeric(df$weight_kgs) else rep(NA_real_, nrow(df))
      height_cm <- if ("height_cm" %in% names(df)) as.numeric(df$height_cm) else rep(NA_real_, nrow(df))
      hc_cm     <- if ("childTimepointHC_MRI_cm" %in% names(df)) as.numeric(df$childTimepointHC_MRI_cm) else rep(NA_real_, nrow(df))
      
      # Initialize columns with NA (They will stay NA if not overwritten)
      for (col in c("zlen", "zwei", "zwfl", "zbmi", "zhc")) {
        if (!col %in% names(df)) df[, (col) := NA_real_]
      }
      
      # --- WHO Transition Gap Fix (60 to 61 months) ---
      gap_idx <- which(calc_age_days > 1826 & calc_age_days < 1857)
      if (length(gap_idx) > 0) {
        calc_age_days[gap_idx] <- 1826
        message("Notice: Temporarily adjusted calculation age to 1826 days (60 months) for ", 
                length(gap_idx), " child(ren) to bridge the WHO 60-61 month gap.")
      }
      
      # --- Subset 1: Under 5 years (<= 1826 calculation days) -> WHO anthro ---
      is_under_5 <- calc_age_days <= 1826 & !is.na(calc_age_days)
      if (any(is_under_5)) {
        u5_idx <- which(is_under_5)
        z_u5 <- suppressWarnings(anthro::anthro_zscores(
          sex = sex_vec[u5_idx],
          age = calc_age_days[u5_idx],
          weight = weight_kg[u5_idx],
          lenhei = height_cm[u5_idx],
          headc = hc_cm[u5_idx] 
        ))
        
        df[u5_idx, zlen := z_u5$zlen]
        df[u5_idx, zwei := z_u5$zwei]
        df[u5_idx, zwfl := z_u5$zwfl]
        df[u5_idx, zbmi := z_u5$zbmi]
        df[u5_idx, zhc := z_u5$zhc] 
      }
      
      # --- Subset 2: Over 5 years (> 1826 calculation days) -> WHO anthroplus ---
      is_over_5 <- calc_age_days > 1826 & !is.na(calc_age_days)
      if (any(is_over_5)) {
        
        if (!requireNamespace("anthroplus", quietly = TRUE)) {
          message("Notice: 'anthroplus' package is required for children > 5 years but is not installed. Skipping older kids.")
        } else {
          o5_idx <- which(is_over_5)
          
          age_months_o5 <- calc_age_days[o5_idx] / 30.4375
          
          z_o5 <- suppressWarnings(anthroplus::anthroplus_zscores(
            sex = sex_vec[o5_idx],
            age_in_months = age_months_o5,
            weight_in_kg = weight_kg[o5_idx],
            height_in_cm = height_cm[o5_idx]
          ))
          
          df[o5_idx, zlen := z_o5$zhfa]
          df[o5_idx, zwei := z_o5$zwfa] 
          df[o5_idx, zbmi := z_o5$zbfa]
        }
      }
      
      # 3. Clean up: Note 100% NA columns but DO NOT drop them!
      skipped <- c()
      if (all(is.na(df$zlen))) { skipped <- c(skipped, "zlen") }
      if (all(is.na(df$zwei))) { skipped <- c(skipped, "zwei") }
      if (all(is.na(df$zwfl))) { skipped <- c(skipped, "zwfl") }
      if (all(is.na(df$zbmi))) { skipped <- c(skipped, "zbmi") }
      if (all(is.na(df$zhc)))  { skipped <- c(skipped, "zhc") } 
      
      if (length(skipped) > 0) {
        message("Notice: The following z-scores evaluated to 100% NA (columns kept as NA to prevent downstream errors): ", paste(skipped, collapse = ", "))
      }
    }
  }
  
  # 4. Define indicators safely
  if ("zlen" %in% names(df)) df[, stunting := fifelse(zlen < -2, 1, 0)]
  if ("zwei" %in% names(df)) df[, underweight := fifelse(zwei < -2, 1, 0)]
  if ("zwfl" %in% names(df)) df[, wasting := fifelse(zwfl < -2, 1, 0)]
  
  # 5. Calculate longitudinal 'ever' indicators safely
  if ("subject_id" %in% names(df)) {
    if ("stunting" %in% names(df)) {
      df[, ever_stunted := if (all(is.na(stunting))) NA_real_ else max(stunting, na.rm = TRUE), by = subject_id]
    }
    if ("underweight" %in% names(df)) {
      df[, ever_underweight := if (all(is.na(underweight))) NA_real_ else max(underweight, na.rm = TRUE), by = subject_id]
    }
    if ("wasting" %in% names(df)) {
      df[, ever_wasted := if (all(is.na(wasting))) NA_real_ else max(wasting, na.rm = TRUE), by = subject_id]
    }
  }
  
  return(df)
}

dt_orig[is.na(childTimepointAge_days)&!is.na(childTimepointAge_months), childTimepointAge_days := childTimepointAge_months * 30.4375]

dt_orig[, weight_kgs := childTimepointWeight_kgs]
dt_orig[, height_cm := childTimepointHeight_inches * 2.54]
#dt_orig[, childBiologicalSex := tools::toTitleCase(tolower(childBiologicalSex))]
dt_orig[, child_sex := childBiologicalSex]
dt_orig[, other_age_days := childTimepointAge_days]
dt_orig[, subject_id := UniqueStudyID]

dt_processed <- process_anthro(dt_orig)

dt_processed[, `:=` (BMI = bmi,
                     BMIZ = zbmi,
                     WAZ = zwei,
                     HAZ = zlen,
                     WHZ = zwfl,
                     HCZ = zhc, 
                     CurrentlyStunted = stunting,
                     CurrentlyWasted = wasting,
                     EverStunted = ever_stunted,
                     EverWasted = ever_wasted)]

dt_processed <- subset(dt_processed, select = -c(weight_kgs, height_cm, child_sex, other_age_days, 
                                                 subject_id, bmi, zbmi, zwei, zlen, zwfl, zhc, stunting, 
                                                 wasting, underweight, ever_stunted, ever_wasted, ever_underweight))

## Calculate maximum age GSED and BSID scores if available ----------------------------------------------------------------------
scores_to_check <- c("gsed_LongForm_DAZScore", "bsid_CCS", "bsid_LCS", "bsid_MCS")

for (score_col in scores_to_check) {
  
  if (score_col %in% names(dt_processed) && "childTimepointAge_days" %in% names(dt_processed)) {
    
    # Define exact output column names
    out_score <- paste0("MAX_", score_col)
    out_age   <- paste0("MAX_", score_col, "_Age_days")
    
    # 1. Force raw data to be numeric (keeps your warnings active!)
    dt_processed[[score_col]] <- as.numeric(dt_processed[[score_col]])
    dt_processed$childTimepointAge_days <- as.numeric(dt_processed$childTimepointAge_days)
    
    # 2. Scrub the existing MAX_ columns of any TRUE/FALSE values by forcing them to numeric
    if (out_score %in% names(dt_processed)) dt_processed[[out_score]] <- as.numeric(dt_processed[[out_score]])
    if (out_age %in% names(dt_processed)) dt_processed[[out_age]] <- as.numeric(dt_processed[[out_age]])
    
    # 3. Create a clean subset to find the maximums
    clean_subset <- dt_processed[!is.na(dt_processed[[score_col]]) & !is.na(dt_processed$childTimepointAge_days), ]
    
    if (nrow(clean_subset) > 0) {
      
      # Sort strictly by ID, then by age
      clean_subset <- clean_subset[order(UniqueStudyID, childTimepointAge_days)]
      
      # Grab the very last row for each ID (which represents their oldest age)
      max_values <- clean_subset[, .SD[.N], by = UniqueStudyID]
      
      # 4. The Magic Bullet: Use match() to slide the data perfectly into the existing columns!
      # This finds the exact row mapping between your main dataset and the max_values table.
      idx <- match(dt_processed$UniqueStudyID, max_values$UniqueStudyID)
      
      # Write the data. If a subject doesn't exist in max_values, it safely leaves an NA.
      dt_processed[[out_score]] <- max_values[[score_col]][idx]
      dt_processed[[out_age]]   <- max_values$childTimepointAge_days[idx]
      
      message("Successfully populated existing max columns for: ", score_col)
      
    } else {
      message("Skipping ", score_col, ": All valid values are NA.")
    }
  }
}

## Save data ----------------------------------------------------------------------
write.csv(dt_processed, paste0(save_dir, file_name), row.names = FALSE)
