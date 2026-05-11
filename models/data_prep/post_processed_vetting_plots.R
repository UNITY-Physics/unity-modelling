##########################################################
## Vetting post-processed cohort data
## Brain Health Metrics 
## April 2026
##########################################################

## Set up ----------------------------------------------------------------------
rm(list=ls())
username <- Sys.getenv("USER")
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(openxlsx)

l_root <- "/FILEPATH/"
read_dir <- paste0(l_root, "FILEPATH/")
save_dir <- paste0(l_root, "FILEPATH/")

extract_country <- "Uganda" 
extract_cohort <- "LONISAC" 


## Read in data and combine ----------------------------------------------------------------------
dt_orig <- as.data.table(read.csv(paste0(read_dir, extract_cohort, "_" , extract_country, "_harmonized_data.csv")))


## Define variables for plotting ----------------------------------------------------------------------

# 1. Variables for Density Plots (Distribution only)
continuous_vars_density <- c("household_size", "household_adults", "household_children", 
                             "mother_educ", "maternal_age_years", "asset_count", "birth_order", "zlen", 
                             "zwei", "zwfl", "zbmi", "gestation_weeks", "bw_grams", 
                             "birthlength_cm", "apgar_1min", "apgar_5min", "epds_score",
                             "gsed_lf_daz", "bsid_ccs", "bsid_lcs", "bsid_mcs",
                             "gsed_lf_daz_max", "bsid_ccs_max", "bsid_lcs_max", "bsid_mcs_max",
                             "maternal_hb", "child_hb")

# 2. Variables for Scatter Plots (Growth/Change over Age)
continuous_vars_overage <- c("height_cm", "weight_kgs", "bmi", "hc_mri_cm",
                             "gsed_lf_daz", "bsid_ccs", "bsid_lcs", "bsid_mcs")

# 3. Categorical Variables (Static)
categorical_vars <- c("child_sex", "num_lang", "electricity", "water", "toilet", 
                      "stove", "telephone", "car", "bicycle", "television", 
                      "floor_material", "wall_material", "roof_material", "stunting", 
                      "underweight", "wasting", "ever_stunted", 
                      "ever_underweight", "ever_wasted", "mother_educ_cat", "preg_high_bp", 
                      "preg_dm", "preg_fever", "delivery_compl", "nicu", "maternal_hiv", 
                      "ever_bf")


## Demographic comparisons ----------------------------------------------------------------------
# Create dynamic output path for the specific cohort
output_path <- paste0(save_dir, extract_cohort, "_", extract_country, "_demographic_vetting.pdf")

# Factor ordering for education if present
if ("mother_educ_cat" %in% names(dt_orig)) {
  dt_orig[, mother_educ_cat := factor(mother_educ_cat, levels = c(
    "No education", 
    "1 to 5 years", 
    "6 to 11 years", 
    "12 years", 
    "Greater than 12 years"
  ))]
}

pdf(output_path, width = 10, height = 7)

## -----------------------------------------------------------------------------
## 1. Density Plots (Continuous Variables)
## -----------------------------------------------------------------------------
present_density_vars <- intersect(names(dt_orig), continuous_vars_density)

for (var in present_density_vars) {
  
  # Filter data
  plot_data <- dt_orig[!is.na(get(var))]
  # Keep only unique subject-level data for static variables
  plot_data <- subset(plot_data, select = c("subject_id", var))
  plot_data <- unique(plot_data)
  
  if (nrow(plot_data) == 0) next
  
  # --- DYNAMIC SS CALCULATION ---
  n_subjects <- uniqueN(plot_data$subject_id)
  ss_caption <- paste0("Subjects with data: N=", n_subjects)
  
  # Calculate Mean
  mean_val <- mean(plot_data[[var]], na.rm = TRUE)
  
  p <- ggplot(plot_data, aes(x = get(var))) +
    geom_histogram(aes(y = after_stat(density)), 
                   fill = "steelblue",
                   alpha = 0.6, 
                   color = "black", 
                   bins = 30) +
    geom_vline(xintercept = mean_val, color = "darkred", 
               linetype = "dashed", size = 1.2) +
    
    labs(title = paste("Distribution of", var),
         subtitle = "Dashed red line indicates mean value",
         caption = ss_caption,
         x = var,
         y = "Density") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          plot.caption = element_text(hjust = 1, face = "italic", color = "grey30"))
  
  print(p)
}

## -----------------------------------------------------------------------------
## 2. Scatter Plots Over Age (Continuous Variables)
## -----------------------------------------------------------------------------
if ("other_age_days" %in% names(dt_orig)) {
  
  present_age_vars <- intersect(names(dt_orig), continuous_vars_overage)
  
  for (var in present_age_vars) {
    
    # Filter data (Must have Variable AND Age)
    plot_data <- dt_orig[!is.na(get(var)) & !is.na(other_age_days)]
    
    if (nrow(plot_data) == 0) next
    
    # --- DYNAMIC SS CALCULATION ---
    n_subjects <- uniqueN(plot_data$subject_id)
    ss_caption <- paste0("Subjects with data: N=", n_subjects)
    
    p <- ggplot(plot_data, aes(x = other_age_days, y = get(var))) +
      geom_point(alpha = 0.5, size = 1.5, color = "steelblue") +
      geom_smooth(method = "loess", color = "darkblue", fill = "lightblue", alpha = 0.3) +
      
      labs(title = paste(var, "vs. Age"),
           subtitle = "Scatter plot with Loess smoothing",
           caption = ss_caption,
           x = "Age (days)",
           y = var) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            plot.caption = element_text(hjust = 1, face = "italic", color = "grey30"))
    
    print(p)
  }
} else {
  message("Skipping scatter plots: 'other_age_days' variable not found in dataset.")
}

## -----------------------------------------------------------------------------
## 3. Bar Plots (Categorical Variables - Static)
## -----------------------------------------------------------------------------
present_cat_vars <- intersect(names(dt_orig), categorical_vars)

for (var in present_cat_vars) {
  
  plot_data <- dt_orig[!is.na(get(var))]
  plot_data <- subset(plot_data, select = c("subject_id", var))
  plot_data <- unique(plot_data)
  
  if (nrow(plot_data) == 0) next
  
  # --- DYNAMIC SS CALCULATION ---
  n_subjects <- uniqueN(plot_data$subject_id)
  ss_caption <- paste0("Subjects with data: N=", n_subjects)
  
  # --- Check for 0/1 Binary Coding and Recode ---
  unique_vals <- unique(plot_data[[var]])
  if (is.numeric(unique_vals) && all(unique_vals %in% c(0, 1), na.rm = TRUE)) {
    plot_data[, (var) := ifelse(get(var) == 1, "Yes", "No")]
  }
  
  # Summarize data
  summary_dt <- plot_data[, .N, by = c(var)]
  setnames(summary_dt, var, "val") # Rename safely for plotting
  summary_dt[, prop := N / sum(N)]
  summary_dt[, label := paste0(round(prop * 100, 1), "%")]
  
  p <- ggplot(summary_dt, aes(x = as.factor(val), y = prop)) +
    geom_bar(stat = "identity", fill = "steelblue", color = "black", show.legend = FALSE) +
    geom_text(aes(label = label), vjust = -0.5, size = 3.5) +
    
    labs(title = paste("Proportion of", var),
         caption = ss_caption,
         x = var,
         y = "Proportion") +
    scale_y_continuous(labels = scales::percent, limits = c(0, max(summary_dt$prop) * 1.15)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.caption = element_text(hjust = 1, face = "italic", color = "grey30"))
  
  print(p)
}

dev.off()
message("PDF saved to: ", output_path)
