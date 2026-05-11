##########################################################
## Standardizing cohort data - helper functions
## Brain Health Metrics IHME
## April 2026
##########################################################

username <- Sys.getenv("USER")

library(anthro, lib.loc = paste0("/FILEPATH/", username, "/"))

pivot_long_lonisac <- function(dt) {
  # Ensure input is a data.table
  dt <- as.data.table(dt)
  
  # 1. Identify column groups --------------------------
  all_names <- names(dt)
  
  # Helper: Matches 'subject_id' OR the specified regex prefix
  get_cols <- function(prefix_pattern) {
    pattern <- paste0("^subject_id$|", prefix_pattern) 
    grep(pattern, all_names, value = TRUE)
  }
  
  # Pass the regex patterns to capture the columns
  # ^m3_ captures anything starting with m3_, including m3_vm3
  cols_t1 <- get_cols("^rec_")
  cols_t2 <- get_cols("^m3_")
  cols_t3 <- get_cols("^m6_")
  cols_t4 <- get_cols("^m12_") # Added 12 mo
  
  # 2. Create subsets ----------------------------------------------------------
  long_t1 <- dt[, ..cols_t1]
  long_t2 <- dt[, ..cols_t2]
  long_t3 <- dt[, ..cols_t3]
  long_t4 <- dt[, ..cols_t4] # Added 12 mo
  
  # 3. Standardize column names ------------------------------------------------
  clean_names <- function(d, prefix_pattern) {
    old <- names(d)
    
    # gsub replaces the matched prefix pattern with nothing ("")
    new <- gsub(prefix_pattern, "", old)
    setnames(d, old, new)
  }
  
  # Remove prefixes exactly as specified:
  clean_names(long_t1, "^rec_")
  clean_names(long_t2, "^m3_vm3|^m3_")
  clean_names(long_t3, "^m6_vm6|^m6_")
  clean_names(long_t4, "^m12_vm12|^m12_") # Added 12 mo
  
  # 4. Add time labels ---------------------------------------------------------
  long_t1[, time := "Enrollment"]
  long_t2[, time := "3 mo"]
  long_t3[, time := "6 mo"]
  long_t4[, time := "12 mo"] # Added 12 mo
  
  # 5. Bind together -----------------------------------------------------------
  long_dt <- rbindlist(
    list(long_t1, long_t2, long_t3, long_t4), # Added long_t4
    fill = TRUE,
    use.names = TRUE
  )
  
  return(long_dt)
}


pivot_long_accra <- function(dt) {
  # Ensure input is a data.table
  dt <- as.data.table(dt)
  
  # 1. Identify column groups --------------------------
  all_names <- names(dt)
  
  # Helper: Matches 'subject_id' OR the specified regex suffix
  get_cols <- function(suffix_pattern) {
    pattern <- paste0("^subject_id$|", suffix_pattern) 
    # ADDED perl = TRUE to allow negative lookbehinds
    grep(pattern, all_names, value = TRUE, perl = TRUE)
  }
  
  # Pass the regex patterns to capture the columns
  # (?<!_) means "not preceded by an underscore"
  # $ signifies the end of the string
  cols_t4  <- get_cols("(?<!_)_4$")
  cols_t6  <- get_cols("(?<!_)_6$")
  cols_t12 <- get_cols("(?<!_)_12$")
  
  # For Enrollment, grab everything that does NOT match those three specific timepoint patterns
  cols_enroll <- grep("(?<!_)(_4|_6|_12)$", all_names, value = TRUE, invert = TRUE, perl = TRUE)
  
  # 2. Create subsets ----------------------------------------------------------
  long_enroll <- dt[, ..cols_enroll]
  long_t4     <- dt[, ..cols_t4]
  long_t6     <- dt[, ..cols_t6]
  long_t12    <- dt[, ..cols_t12]
  
  # 3. Standardize column names ------------------------------------------------
  clean_names <- function(d, suffix_pattern) {
    old <- names(d)
    
    # gsub replaces the matched suffix pattern with nothing ("")
    # ADDED perl = TRUE here as well
    new <- gsub(suffix_pattern, "", old, perl = TRUE)
    setnames(d, old, new)
  }
  
  # Remove suffixes exactly as specified:
  clean_names(long_t4, "(?<!_)_4$")
  clean_names(long_t6, "(?<!_)_6$")
  clean_names(long_t12, "(?<!_)_12$")
  
  # 4. Add time labels ---------------------------------------------------------
  long_enroll[, time := "Enrollment"]
  long_t4[, time := "3-4 mo"]
  long_t6[, time := "6 mo"]
  long_t12[, time := "12 mo"]
  
  # 5. Bind together -----------------------------------------------------------
  long_dt <- rbindlist(
    list(long_enroll, long_t4, long_t6, long_t12), 
    fill = TRUE,
    use.names = TRUE
  )
  
  return(long_dt)
}


standardize_cohort_data <- function(df, target_filepath, mapping_table) {
  
  # 1. Ensure inputs are data.tables
  df <- as.data.table(df)
  mapping_table <- as.data.table(mapping_table)
  
  # 2. Subset mapping
  target_path_str <- as.character(target_filepath)
  
  cohort_map <- mapping_table[as.character(cohort_filepath) == target_path_str, ]
  
  if (nrow(cohort_map) == 0) {
    stop(paste("No mapping found for:", cohort_filepath))
  }
  
  # --- Identify and Message Missing Variables ---
  expected_cols <- cohort_map$original_name
  actual_cols <- names(df)
  missing_vars <- setdiff(expected_cols, actual_cols)
  
  if (length(missing_vars) > 0) {
    message("-------------------------------------------------------")
    message(paste("NOTE:", length(missing_vars), "variables from the mapping sheet were NOT found in the dataset:"))
    message(paste(missing_vars, collapse = ", "))
    message("-------------------------------------------------------")
  } else {
    message("Success: All mapped variables were found in the dataset.")
  }
  
  # Continue with valid map (intersection only)
  valid_map <- cohort_map[original_name %in% names(df)]
  
  if (nrow(valid_map) == 0) {
    stop("None of the mapped columns exist in the provided dataframe.")
  }
  
  # 3. Subset columns (INCLUDES time)
  # ---------------------------------------------------------
  cols_to_keep <- valid_map$original_name
  
  # Check for 'time' and add to extraction list if present
  has_timepoint <- "time" %in% names(df)
  if (has_timepoint) {
    cols_to_keep <- unique(c(cols_to_keep, "time"))
  }
  
  # Extract columns
  df_subset <- df[, ..cols_to_keep]
  
  # 4. Rename columns (EXCLUDES time)
  map_old <- valid_map$original_name
  map_new <- as.character(valid_map$standard_name)
  
  setnames(df_subset, old = map_old, new = map_new)
  
  # 5. Enforce variable types
  for (i in 1:nrow(valid_map)) {
    col <- map_new[i]
    type <- valid_map$type[i]
    
    if (col %in% names(df_subset)) {
      if (type == "numeric") {
        # FIX: Wrapped get(col) in as.character() to safely handle factors
        df_subset[, (col) := suppressWarnings(as.numeric(as.character(get(col))))]
      } else if (type == "character") {
        df_subset[, (col) := as.character(get(col))]
      } else if (type == "integer") {
        # FIX: Wrapped get(col) in as.character() to safely handle factors
        df_subset[, (col) := suppressWarnings(as.integer(as.character(get(col))))]
      }
    }
  }
  
  return(df_subset)
}


convert_missing_units <- function(df) {
  # Ensure input is a data.table
  df <- as.data.table(df)
  
  # --- Helper Function ---
  convert_if_needed <- function(d, target_col, source_col, factor) {
    
    # Only proceed if source column exists
    if (source_col %in% names(d)) {
      
      # Convert source to numeric (suppress warnings for non-numeric text)
      # FIX: Wrapped in as.character() to safely handle factors
      vals <- suppressWarnings(as.numeric(as.character(d[[source_col]])))
      
      # Scenario A: Target doesn't exist -> Create it (Length N -> Length N)
      if (!target_col %in% names(d)) {
        d[, (target_col) := vals * factor]
      } 
      # Scenario B: Target exists -> Fill missing gaps (NAs) only
      else {
        d[is.na(get(target_col)) & !is.na(vals), 
          (target_col) := vals[.I] * factor]
      }
    }
    return(d)
  }
  
  # --- Wrapper for Bi-Directional Conversion ---
  convert_pair <- function(d, unit_a, unit_b, factor_a_to_b) {
    # Convert A -> B
    d <- convert_if_needed(d, unit_b, unit_a, factor_a_to_b)
    # Convert B -> A (using 1/factor)
    d <- convert_if_needed(d, unit_a, unit_b, 1/factor_a_to_b)
    return(d)
  }
  
  # --- 1. Age Conversions ---
  # Days <-> Months (30.4375 days per month)
  df <- convert_pair(df, "other_age_months", "other_age_days", 30.4375)
  
  # --- 2. Anthropometry ---
  # Height: Inches * 2.54 = cm
  df <- convert_pair(df, "height_inches", "height_cm", 2.54)
  
  # Weight: Kgs * 2.20462 = Lbs
  df <- convert_pair(df, "weight_kgs", "weight_lbs", 2.20462)
  
  # Head Circumference: Inches * 2.54 = cm
  df <- convert_pair(df, "hc_measured_inches", "hc_measured_cm", 2.54)
  
  # --- 3. Gestation ---
  # Weeks * 7 = Days
  df <- convert_pair(df, "gestation_weeks", "gestation_days", 7)
  
  # --- 4. Birth Weight ---
  # Step A: Kgs * 1000 = Grams
  df <- convert_pair(df, "bw_kgs", "bw_grams", 1000)
  
  # Step B: Kgs * 2.20462 = Lbs
  df <- convert_pair(df, "bw_kgs", "bw_lbs", 2.20462)
  
  # Step C: Re-run Grams check (in case Lbs -> Kgs happened, now need Kgs -> Grams)
  df <- convert_pair(df, "bw_kgs", "bw_grams", 1000)
  
  # --- 5. Birth Length ---
  df <- convert_pair(df, "birthlength_inches", "birthlength_cm", 2.54)
  
  # --- 6. MRI HC ---
  df <- convert_pair(df, "hc_mri_inches", "hc_mri_cm", 2.54)
  
  return(df)
}


process_anthro <- function(df) {
  
  # 0. Calculate raw BMI directly safely
  if (all(c("weight_kgs", "height_cm") %in% names(df))) {
    temp_bmi <- as.numeric(df$weight_kgs) / ((as.numeric(df$height_cm) / 100)^2)
    
    # Only add the column if it's not 100% NA
    if (!all(is.na(temp_bmi))) {
      df[, bmi := temp_bmi]
    } else {
      message("Notice: Calculated 'bmi' would be entirely NA. Column not added.")
    }
  } else {
    message("Notice: 'weight_kgs' or 'height_cm' missing. Skipping raw BMI calculation.")
  }
  
  # 1. Check for absolute minimum required variables for WHO z-scores
  if (!all(c("child_sex", "other_age_days") %in% names(df))) {
    message("Notice: Missing 'child_sex' or 'other_age_days'. Skipping WHO z-score calculations entirely.")
    return(df)
  }
  
  # Only run if output columns don't already exist
  if (!all(c("zlen", "zwei", "zwfl", "zbmi") %in% names(df))) {
    
    # Extract core vectors
    sex_vec <- fcase(
      df$child_sex == "Male", 1L,
      df$child_sex == "Female", 2L,
      default = NA_integer_
    )
    age_days <- as.numeric(df$other_age_days)
    
    if (all(is.na(sex_vec)) || all(is.na(age_days))) {
      message("Notice: All rows have missing sex or age data. Skipping WHO z-score calculations.")
    } else {
      
      # Safely pull weight and height if they exist, otherwise fill with NA
      weight_kg <- if ("weight_kgs" %in% names(df)) as.numeric(df$weight_kgs) else rep(NA_real_, nrow(df))
      height_cm <- if ("height_cm" %in% names(df)) as.numeric(df$height_cm) else rep(NA_real_, nrow(df))
      
      # 2. Calculate Z-scores (suppress warnings in case we feed it 100% NAs for weight/height)
      zscores <- suppressWarnings(anthro_zscores(
        sex = sex_vec,
        age = age_days,
        weight = weight_kg,
        lenhei = height_cm
      ))
      zscores <- as.data.table(zscores)
      
      # 3. Bind results back ONLY if they contain at least one valid number
      skipped <- c()
      
      if (!all(is.na(zscores$zlen))) df[, zlen := zscores$zlen] else skipped <- c(skipped, "zlen")
      if (!all(is.na(zscores$zwei))) df[, zwei := zscores$zwei] else skipped <- c(skipped, "zwei")
      if (!all(is.na(zscores$zwfl))) df[, zwfl := zscores$zwfl] else skipped <- c(skipped, "zwfl")
      if (!all(is.na(zscores$zbmi))) df[, zbmi := zscores$zbmi] else skipped <- c(skipped, "zbmi")
      
      # Let the user know if any scores evaluated to 100% NA and were dropped
      if (length(skipped) > 0) {
        message("Notice: The following z-scores evaluated to 100% NA and were not added: ", paste(skipped, collapse = ", "))
      }
    }
  }
  
  # 4. Define indicators safely
  # Because these check if the z-score column exists, they will automatically 
  # skip creation if the z-scores were dropped in the step above!
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


