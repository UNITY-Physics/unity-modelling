##########################################################
## Standardizing cohort data for analyses - for FlyWheel upload
## Brain Health Metrics IHME
## April 2026
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

#System configuration
l_root <- "/FILEPATH/"
j_root <- "/FILEPATH/"

#Set objects
mapping <- as.data.table(read.csv("/FILEPATH/cohort_data_mapping.csv"))
var_map <- fread(paste0("/FILEPATH/variable_mapping.csv"))

extract_country <- "Uganda" # UPDATE
extract_cohort <- "LONISAC" # UPDATE

mapping <- mapping[country == extract_country & cohort_name == extract_cohort, ]

filename <- mapping$cohort_filename

read_dir <- mapping$cohort_filepath[1]
save_dir <- paste0(l_root, "FILEPATH/")

code_dir <- paste0("/FILEPATH/", username, "/early_brain_health_cohorts/")
source(paste0(code_dir, "process_cohort_data_unity/standardize_cohort_data_functions.R"))


## Read cohort data ----------------------------------------------------------------------
# 1. Extract the column name to merge by if multiple different datasets
merge_cols <- unique(mapping$merge_column)

# 2. Read in all datasets into a list
df_list <- list()

for (i in 1:nrow(mapping)) {
  # Get row-specific file details
  file_n <- mapping[i, cohort_filename]
  dir_n  <- mapping[i, cohort_filepath]
  full_path <- paste0(dir_n, file_n)
  
  # Read in data with case-insensitive extension matching
  if (grepl("\\.csv$", file_n, ignore.case = TRUE)) {
    temp_df <- fread(full_path, encoding = "Latin-1")
  } else if (grepl("\\.xlsx?$", file_n, ignore.case = TRUE)) {
    temp_df <- as.data.table(read_excel(full_path))
  } else if (grepl("\\.dta$", file_n, ignore.case = TRUE)) {
    temp_df <- haven::read_dta(full_path)
    temp_df <- haven::zap_labels(temp_df)
    temp_df <- haven::zap_label(temp_df)
    temp_df <- haven::zap_formats(temp_df)
    temp_df <- as.data.table(temp_df)
  } else {
    stop(paste("Unsupported file extension for:", file_n))
  }
  
  # --- Standardize the merge column name ---
  current_merge_col <- mapping[i, merge_column]
  
  if (!current_merge_col == "") {
    if (current_merge_col %in% names(temp_df)) {
      # Rename whatever the local ID column is to a universal "subject_id"
      setnames(temp_df, current_merge_col, "subject_id")
    } else {
      stop(paste("Error: Merge column '", current_merge_col, "' not found in file", file_n))
    }
  }
  
  # --- Add custom prefix to all other columns ---
  current_prefix <- mapping[i, add_prefix]
  
  # Check if a valid prefix exists for this specific file (not NA and not empty)
  if (!is.na(current_prefix) && trimws(current_prefix) != "") {
    
    # Grab all column names EXCEPT the merge column
    cols_to_prefix <- setdiff(names(temp_df), "subject_id")
    
    # Create the new names
    new_colnames <- paste0(current_prefix, cols_to_prefix)
    
    # Apply them back to the data.table
    setnames(temp_df, old = cols_to_prefix, new = new_colnames)
  }
  
  # --- Save individual dataset to the global environment ---
  # This creates objects named df_1, df_2, etc., in your workspace
  obj_name <- paste0("df_", i)
  assign(obj_name, temp_df, envir = .GlobalEnv)
  
  message("Saved separate object: ", obj_name, " (", file_n, ")")
  
  # Add to the list for merging
  df_list[[i]] <- temp_df
}

# 3. Merge the datasets together
if (length(df_list) == 1) {
  # Only one dataset; no merging needed
  df_orig <- df_list[[1]]
} else {
  # Multiple datasets; merge them iteratively
  
  # Reduce iteratively applies the merge function across all items in the list
  # We now merge specifically by our newly standardized "subject_id" column
  df_orig <- Reduce(function(x, y) merge(x, y, by = "subject_id", all = TRUE), df_list)
}

if (extract_cohort == "LONISAC" & extract_country == "Uganda") {
  df_orig <- df_orig[!is.na(subject_id),] #some blank rows
} 

if (extract_cohort == "BRAC" & extract_country == "Bangladesh") {
  df_orig <- df_orig[!is.na(session_id),] #only subjects with MRI scan
} 

if (extract_cohort == "REVAMP" & extract_country == "Malawi") {
  for (col in names(df_orig)) {
    # Find the row indices where the column equals "."
    rows_to_replace <- which(df_orig[[col]] == ".")
    
    # If any exist, replace them with NA in place
    if (length(rows_to_replace) > 0) {
      set(df_orig, i = rows_to_replace, j = col, value = NA)
    }
  }
}


## Cleaning/standardizing ----------------------------------------------------------------------
#Convert long if needed using custom functions
if (extract_cohort %in% c("LONISAC") & extract_country == "Uganda") {
  df_orig <- pivot_long_lonisac(df_orig)
}

if (extract_cohort %in% c("Accra") & extract_country == "Ghana") {
  df_orig <- pivot_long_accra(df_orig)
}

#Fix skip patterns and other coding 
if (extract_cohort == "LONISAC" & extract_country == "Uganda" ) {
  df_orig[b01g06 %in% c("1", "Male"), b01g06 := "Male"]
  df_orig[b01g06 %in% c("2", "Female"), b01g06 := "Female"]
  
  df_orig[, a00f02 := fcase(a00f02 == "Yes", 1, a00f02 == "No", 0, default = NA_real_)]
  
  df_orig[, d01 := fcase(d01 %in% c("1", "Yes"), 1, d01 %in% c("0", "No"), 0, default = NA_real_)]
  df_orig[, d02 := fcase(d02 %in% c("1", "Yes"), 1, d02 %in% c("0", "No"), 0, default = NA_real_)]
  df_orig[d01 == 1, d02 := 1]
  
  df_orig[subject_id == 1163 & time == "Enrollment", `:=` (b01g02a = 50, b01g02b = 50)]
  df_orig[subject_id == 1198 & time == "Enrollment", `:=` (b01g02a = 52, b01g02b = 52)]
  df_orig[subject_id == 1081 & time == "Enrollment", `:=` (b01g02a = 47.3, b01g02b = 47.3)]
  
  df_orig[, f02 := (f02a+f02b)/2]
  df_orig[, f03 := (f03a+f03b)/2]
  df_orig[time == "Enrollment", b01g03 := (b01g03a+b01g03b)/2]
  df_orig[time == "Enrollment", b01g02 := (b01g02a+b01g02b)/2]
  
  df_orig[a00f04 %in% c("1"), a00f04 := "Wood"]
  df_orig[a00f04 %in% c("2"), a00f04 := "Charcoal"]
  df_orig[a00f04 %in% c("4"), a00f04 := "Gas"]

  df_orig[a00h04 %in% c("1", "Earth/Dung"), a00h04 := "Earth/dung/sand"]
  df_orig[a00h04 %in% c("2", "Cement"), a00h04 := "Cement"]
  df_orig[a00h04 %in% c("3", "Tiles"), a00h04 := "Tiles"]
  df_orig[a00h04 %in% c("6", "Carpet /Vinyl"), a00h04 := "Other"]
  df_orig[a00h04 %in% c("7", "Other, specify_________"), a00h04 := "Other"]
  
  df_orig[a00h05 %in% c("1", "Thatch, grass"), a00h05 := "Thatch/grass"]
  df_orig[a00h05 %in% c("2", "Iron sheets"), a00h05 := "Iron sheets"]
  df_orig[a00h05 %in% c("3", "Tiles"), a00h05 := "Tiles"]
  df_orig[a00h05 %in% c("Other, specify______"), a00h05 := "Other"]
  
  df_orig[a00h06 %in% c("4", "5", "Bricks without mortar", "Burnt brick with mortar"), a00h06 := "Brick"]
  df_orig[a00h06 %in% c("1", "Mud and pole"), a00h06 := "Mud"]
  df_orig[a00h06 %in% c("Wood"), a00h06 := "Wood"]
  df_orig[a00h06 %in% c("6", "Plastered walls"), a00h06 := "Plastered walls"]
  df_orig[a00h06 %in% c("7", "Other, specify___"), a00h06 := "Other"]
  
  df_orig[a00h02 %in% c("1", "Open pit"), a00h02 := "Open pit"]
  df_orig[a00h02 %in% c("2", "Pit latrines"), a00h02 := "Pit latrines"]
  df_orig[a00h02 %in% c("3", "VIP latrine"), a00h02 := "VIP latrine"]
  df_orig[a00h02 %in% c("4", "Flush toilet"), a00h02 := "Flush toilet"]
  df_orig[a00h02 %in% c("Other, specify__"), a00h02 := "Other"]
  df_orig[a00h02 %in% c("0"), a00h02 := NA] # no coding for 0
  
  df_orig[, a04 := fcase(a04 %in% c("1", "Yes"), 1, a04 %in% c("0", "No"), 0, default = NA_real_)]
  
  df_orig[, water := fifelse(is.na(a00f05), NA_real_, 
                             fifelse(a00f05 %like% "Piped into dwelling" | a00f05 %like% "Piped into yard / plot", 1, 0))]
  df_orig[, toilet := fifelse(is.na(a00h02), NA_real_, 
                             fifelse(a00h02 %like% "Flush toilet", 1, 0))]
  
  df_orig[b01o02>=1, b01o02 := 1]
  
  df_orig[epds_score>100, epds_score:=NA] #impossible scores
  
  df_orig[!time == "Enrollment", a00c01 := NA]
  df_orig[a00c01%in%c(15:18), a00c01:=13+a00c01_years] #adjusting for years in tertiary education
  df_orig[a00c01%in%c(19), a00c01:=NA] #"don't know" response for maternal education
  
  df_orig[time == "6 mo", maternal_HB := NA] #not required for 6 mo visit so lots of missingness; non-missing values are rounded to whole number
  
}

if (extract_cohort == "Accra" & extract_country == "Ghana" ) {
  df_orig[, baby_sex := as.character(baby_sex)]
  df_orig[baby_sex == "0", baby_sex := "Male"]
  df_orig[baby_sex == "1", baby_sex := "Female"]
  
  df_orig[, excl_bf_6mo := {
    # Extract the mode_feed value specifically at the 6-month timepoint
    val <- mode_feed[time == "6 mo"]
    
    # Determine the true/false/NA status based on the 6 mo data
    bf_val <- if (length(val) == 0 || is.na(val[1])) {
      NA_real_
    } else if (val[1] == 1) {
      1
    } else {
      0
    }
    
    # Return the value ONLY for the Enrollment row, otherwise leave as NA
    fifelse(time == "Enrollment", bf_val, NA_real_)
    
  }, by = subject_id]
  
  df_orig[, mother_educ := fcase( 
    educ == 0, 0,   # None
    educ == 1, 6,   # Primary
    educ == 2, 9,   # JHS 
    educ == 3, 12,  # SHS 
    educ == 4, 16,  # Tertiary
    default = NA_real_
  )]
  
  #calculating age at timepoint using dob and date of timepoint
  # 1. Convert both columns to formal Date objects
  df_orig[, dob := as.Date(dob, format = "%Y-%m-%d")]
  df_orig[, date := as.Date(date, format = "%Y-%m-%d")]
  
  # 2. Broadcast the Enrollment 'dob' to all rows for each subject
  # This groups by subject_id, finds the non-NA dob, and fills it everywhere
  df_orig[, dob := dob[!is.na(dob)][1], by = subject_id]
  
  # 3. Calculate other_age_days for the specific timepoints
  df_orig[time %in% c("3-4 mo", "6 mo", "12 mo") & !is.na(dob) & !is.na(date), 
          other_age_days := as.numeric(date - dob)]
  
  #fixing some birthweight values that are grams instead of kg
  df_orig[baby_weight > 100, baby_weight := baby_weight / 1000]
  
  #anthro error fixes
  df_orig[subject_id == "HYPERFINE-0067"&time=="3-4 mo"&lenght==545, lenght := 54.5] #fixing a length value - assuming missing a decimal
  df_orig[subject_id == "HYPERFINE-0108"&time=="3-4 mo"&lenght==46, lenght := NA] #value implausibly low and is equal to enrollment length, so setting to missing
  df_orig[subject_id == "HYPERFINE-0116"&time=="6 mo"&lenght==47, lenght := NA] #value implausibly low and is less than enrollment and 3-4 mo length, so setting to missing
  
  #using weight and length at enrollment as birth weight and length
  #df_orig[time == "Enrollment", `:=` (bw_kgs = weight, birthlength_cm = lenght)]
  
}

if (extract_cohort == "BRAC" & extract_country == "Bangladesh" ) {
  df_orig[, maternal_age_years := round(b6_1-(childage/12))]
  
  df_orig[, b11_1 := as.numeric(b11_1)]
  df_orig[b11_1 == 51, mother_educ := 0] # did not pass a class
  df_orig[b11_1 == 10, mother_educ := 10] #SSC / Dakhil 
  df_orig[b11_1 == 11, mother_educ := 12] #HSC / Alim  
  df_orig[b11_1 == 12, mother_educ := 14] #B.A / B.Sc / B.Com / Fazil 
  df_orig[b11_1 == 13, mother_educ := 16] #M.A / M.Sc / M.Com / Kamil 
  df_orig[b11_1 == 14, mother_educ := 14] #Diploma / Vocational 
  df_orig[b11_1 == 15, mother_educ := 5] #Hafiz
  df_orig[b11_1 == 50, mother_educ := 3] #Religious education
  df_orig[!is.na(b11_1)&is.na(mother_educ), mother_educ := b11_1] #for values 0 to 9, which represent years of education, just carry those forward

  df_orig[, child_gender := as.character(child_gender)]
  df_orig[child_gender == "1", child_gender := "Male"]
  df_orig[child_gender == "2", child_gender := "Female"]
  
}

if (extract_cohort == "MINE_HF" & extract_country == "Pakistan" ) {
  df_orig[, water := fifelse(household_water_source == "Piped water into dwelling", 1, 0)]
  
  df_orig[, telephone := fcase(
    family_owns_cellphone == "Yes", 1, 
    family_owns_cellphone == "No", 0, 
    default = NA_real_
  )]
  
  df_orig[EPDS_Score == 99, EPDS_Score := NA] #impossible score
  
  # df_orig[maternalEducation_schoolingYears == "0 years", mother_educ := 0]
  # df_orig[maternalEducation_schoolingYears == "0–3 years", mother_educ := 3]
  # df_orig[maternalEducation_schoolingYears == "3–11 years", mother_educ := 11]
  # df_orig[maternalEducation_schoolingYears == "11–13 years", mother_educ := 13]
  # df_orig[maternalEducation_schoolingYears == ">15 years", mother_educ := 15]
  # df_orig[maternalEducation_schoolingYears == "Others", mother_educ := NA]
  
  df_orig[maternalEducation_schoolingYears == "0 years", mother_educ := 0]
  df_orig[maternalEducation_schoolingYears == "0–3 years", mother_educ := 1.5]
  df_orig[maternalEducation_schoolingYears == "3–11 years", mother_educ := 7]
  df_orig[maternalEducation_schoolingYears == "11–13 years", mother_educ := 12]
  df_orig[maternalEducation_schoolingYears == ">15 years", mother_educ := 15]
  df_orig[maternalEducation_schoolingYears == "Others", mother_educ := NA]
  
  df_orig[, household_children := No_children+No_Adolescents]
  df_orig[, household_adults := Female19.x+Male19.x]
  
}

if (extract_cohort == "REVAMP" & extract_country == "Malawi" ) {
  df_orig[maternal_age == "<=16", maternal_age := 16]
  df_orig[maternal_age == ">=35", maternal_age := 35]
  df_orig[, maternal_age := as.numeric(maternal_age)]
  
  df_orig[, age_at_scan := gsub(" months", " mo", age_at_scan)]
  
  df_orig[age_at_scan == "3 mo", infant_weight := as.numeric(infant_weight_3mths)/1000]
  df_orig[age_at_scan == "12 mo", infant_weight := as.numeric(infant_weight_12mths)/1000]
  
  df_orig[age_at_scan == "3 mo", infant_height := as.numeric(infant_height_3mths)]
  df_orig[age_at_scan == "12 mo", infant_height := as.numeric(infant_height_12mths)]
  
  df_orig[age_at_scan == "3 mo", infant_hc := as.numeric(infant_hc_3mths)]
  df_orig[age_at_scan == "12 mo", infant_hc := as.numeric(infant_hc_12mths)]
  
  df_orig[age_at_scan == "12 mo", BSID_cognition := as.numeric(BSID_cognition_12mths)]
  df_orig[age_at_scan == "12 mo", BSID_language := as.numeric(BSID_language_12mths)]
  df_orig[age_at_scan == "12 mo", BSID_motor := as.numeric(BSID_motor_12mths)]
  
}

#Standardize and select column names
df_clean <- standardize_cohort_data(df = df_orig, 
                                    target_filepath = read_dir,
                                    mapping_table = var_map)


# Separate and merge one-time measured variables if needed
if ((extract_cohort %in% c("LONISAC") & extract_country == "Uganda")|
    (extract_cohort %in% c("Accra") & extract_country == "Ghana")) {
  
  # 1. Define One-Time Variables (Potential List)
  potential_onetime_vars <- c("subject_id", 
                              "child_sex", "household_size", "household_adults", "household_children", 
                              "gestation_weeks", "bw_grams", "bw_kgs", "birthlength_cm", "maternal_age_years",
                              "maternal_hiv", "electricity", "water", "toilet", "stove", 
                              "telephone", "car", "bicycle", "television", "floor_material", "wall_material", 
                              "roof_material", "main_cooking_fuel", "main_water_source", "toilet_type", "num_lang", 
                              "trimester1_alc_binary", "trimester2_alc_binary", "trimester3_alc_binary", 
                              "preg_high_bp", "preg_dm", "preg_fever", "num_pregnancies", "anc_visits","anc_start_gestation_weeks",
                              "apgar_1min", "apgar_5min", "delivery_compl", "nicu", "jaundice", 
                              "ever_bf", "preg_min_hb", "preg_min_hb_ga", "preg_min_hb_trimester", 
                              "preg_anemia_status", "preg_anemia_severity", "mother_educ", "excl_bf_6mo", "birth_order")
  
  # Ensure subject_id is always kept
  present_vars <- intersect(potential_onetime_vars, names(df_clean))
  
  # Proceed only if we have variables to process (beyond just subject_id)
  if (length(present_vars) > 1) {
    
    if ("time" %in% names(df_clean)) {
      
      # 2a. Extract Baseline (Enrollment) data for these specific variables
      df_baseline_raw <- df_clean[time == "Enrollment", ..present_vars]
      
      # 2b. Consolidate baseline in case of accidental duplicates (safeguard)
      df_baseline <- df_baseline_raw[, lapply(.SD, function(x) {
        val <- x[!is.na(x)]
        if (length(val) > 0) return(val[1]) else return(x[1]) 
      }), by = subject_id]
      
      # Rename baseline columns so they don't clash during the merge
      cols_to_fill <- setdiff(present_vars, "subject_id")
      setnames(df_baseline, old = cols_to_fill, new = paste0("base_", cols_to_fill))
      
      # 3. Remove Enrollment rows from the main dataset
      df_clean <- df_clean[time != "Enrollment", ]
      
      # 4. CRITICAL FIX: Remove empty "ghost" rows BEFORE merging!
      # We check the remaining 3mo and 6mo rows. If everything except ID and time is NA, delete the row.
      check_cols <- setdiff(names(df_clean), c("subject_id", "time"))
      
      if (length(check_cols) > 0) {
        df_clean <- df_clean[rowSums(!is.na(df_clean[, ..check_cols])) > 0, ]
      }
      
      # 5. Merge baseline values onto the VALID non-enrollment rows
      df_clean <- merge(df_clean, df_baseline, by = "subject_id", all.x = TRUE)
      
      # 6. COALESCE: Fill NAs in the current timepoint with the baseline value
      for (col in cols_to_fill) {
        base_col <- paste0("base_", col)
        
        # If the column exists longitudinally, fill the NAs. If it doesn't, create it.
        if (col %in% names(df_clean)) {
          df_clean[is.na(get(col)), (col) := get(base_col)]
        } else {
          df_clean[, (col) := get(base_col)]
        }
        
        # Drop the temporary base_ column
        df_clean[, (base_col) := NULL]
      }
      
      message("Verification passed: Missing longitudinal data filled with Enrollment baseline values.")
      
    } else {
      message("Warning: 'time' variable not found. Cannot perform baseline carry-forward.")
    }
    
  } else {
    message("Warning: None of the specified one-time variables were found in the dataset. Skipping baseline fill step.")
  }
}

#Unit conversions
df_clean <- convert_missing_units(df_clean)

#Process anthropometrics
df_clean <- process_anthro(df_clean)

#Add maternal edu categorical variable
if ("mother_educ" %in% names(df_clean)) {
df_clean[, mother_educ_cat := fifelse(mother_educ == 0, "No education",
                                      fifelse(mother_educ >0 & mother_educ < 6, "1 to 5 years",
                                              fifelse(mother_educ >=6 & mother_educ < 12, "6 to 11 years",
                                                      fifelse(mother_educ == 12, "12 years",
                                                              fifelse(mother_educ > 12, "Greater than 12 years", NA_character_)))))]
}

#Add maternal and child anemia status variables if not present
if ("child_hb" %in% names(df_clean) & !"child_anemia_status" %in% names(df_clean)) {
  df_clean[other_age_months<24&child_hb <10.5, child_anemia_status := 1]
  df_clean[other_age_months<24&child_hb >=10.5, child_anemia_status := 0]
  
  df_clean[other_age_months>=24&child_hb <11, child_anemia_status := 1]
  df_clean[other_age_months>=24&child_hb >=11, child_anemia_status := 0]
}

if ("maternal_hb" %in% names(df_clean) & !"maternal_anemia_status" %in% names(df_clean)) {
  df_clean[maternal_hb <12, maternal_anemia_status := 1]
  df_clean[maternal_hb >=12, maternal_anemia_status := 0]
}

#Add anemia trimester in pregnancy variable
if ("preg_anemia_status" %in% names(df_clean) & "preg_min_hb_trimester" %in% names(df_clean)) {
df_clean[preg_anemia_status == 1, preg_anemia_trimester := preg_min_hb_trimester]
}

#Fix anemia severity coding
if ("preg_anemia_status" %in% names(df_clean) & "preg_anemia_severity" %in% names(df_clean)) {
  df_clean[preg_anemia_status == 0, preg_anemia_severity := 0]
}

#Asset count
if (all(c("electricity", "water", "toilet", "stove", "telephone", "car", "bicycle") %in% names(df_clean))) {
  df_clean[, asset_count := electricity + water + toilet + stove + telephone + car + bicycle]
}

#Calculate maximum age GSED and BSID scores if available
# Define the raw longitudinal score names and their target "max" names
score_mapping <- list(
  "gsed_lf_daz" = "gsed_lf_daz_max",
  "bsid_ccs" = "bsid_ccs_max",
  "bsid_lcs" = "bsid_lcs_max",
  "bsid_mcs" = "bsid_mcs_max"
)

for (raw_score in names(score_mapping)) {
  max_score_name <- score_mapping[[raw_score]]
  max_age_name <- paste0(max_score_name, "_age_days")
  
  # Check if raw data exists AND the max score hasn't already been created
  if (all(c(raw_score, "other_age_days") %in% names(df_clean)) && !(max_score_name %in% names(df_clean))) {
    
    # Subset to rows where neither the score nor the age is missing
    valid_dt <- df_clean[!is.na(get(raw_score)) & !is.na(other_age_days)]
    
    if (nrow(valid_dt) > 0) {
      # Find the row index of the maximum age for each subject
      max_idx <- valid_dt[, .I[which.max(other_age_days)], by = subject_id]$V1
      
      # Extract just the ID, score, and age at that max index
      max_dt <- valid_dt[max_idx, .(subject_id, get(raw_score), other_age_days)]
      
      # Rename the columns to match your standard naming convention
      setnames(max_dt, c("V2", "other_age_days"), c(max_score_name, max_age_name))
      
      # Merge the new static max variables back onto the main dataset
      df_clean <- merge(df_clean, max_dt, by = "subject_id", all.x = TRUE)
      
      message("Successfully calculated terminal maximums for: ", raw_score)
    }
  }
}

# --- Biological Plausibility Checks ---
message("Running biological plausibility checks...")

# 1. Check if birth weight is greater than or equal to follow-up weight
if (all(c("bw_kgs", "weight_kgs", "time") %in% names(df_clean))) {
  bw_violations <- df_clean[time != "Enrollment" & bw_kgs >= weight_kgs, 
                            .(subject_id, time, bw_kgs, weight_kgs)]
  
  if (nrow(bw_violations) > 0) {
    warning(paste("DATA INTEGRITY WARNING: Found", nrow(bw_violations), 
                  "follow-up rows where bw_kgs is >= weight_kgs."))
    print(head(bw_violations)) # Prints the first few offending rows
  } else {
    message("  - Birth weight check passed.")
  }
}

# 2. Check if birth length is greater than or equal to follow-up length
if (all(c("birthlength_cm", "height_cm", "time") %in% names(df_clean))) {
  bl_violations <- df_clean[time != "Enrollment" & birthlength_cm >= height_cm, 
                            .(subject_id, time, birthlength_cm, height_cm)]
  
  if (nrow(bl_violations) > 0) {
    warning(paste("DATA INTEGRITY WARNING: Found", nrow(bl_violations), 
                  "follow-up rows where birthlength_cm is >= height_cm."))
    print(head(bl_violations))
  } else {
    message("  - Birth length check passed.")
  }
}

# 3. Check if child_sex is consistent across timepoints for each subject
if (all(c("subject_id", "child_sex", "time") %in% names(df_clean))) {
  
  # Find subjects with more than one unique non-NA sex value
  sex_violations_ids <- df_clean[!is.na(child_sex), 
                                 .(n_unique = uniqueN(child_sex)), 
                                 by = subject_id][n_unique > 1, subject_id]
  
  if (length(sex_violations_ids) > 0) {
    warning(paste("DATA INTEGRITY WARNING: Found", length(sex_violations_ids), 
                  "subject(s) with conflicting child_sex values across timepoints."))
    
    # Pull the specific rows to show the user exactly where the mismatch is
    sex_violations <- df_clean[subject_id %in% sex_violations_ids, 
                               .(subject_id, time, child_sex)][order(subject_id, time)]
    
    # REMOVED head() so it prints the entire list of violations
    print(sex_violations) 
  } else {
    message("  - Child sex consistency check passed.")
  }
}


## Save Output -------------------------------------------
write.csv(df_clean, paste0(save_dir, extract_cohort, "_" , extract_country, "_harmonized_data.csv"), row.names = FALSE)

