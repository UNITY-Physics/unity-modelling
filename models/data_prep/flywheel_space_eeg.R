######################################################################
## Create SPACE EEG dataset for FlyWheel
## Brain Health Metrics 
## IHME
## August 2025
#####################################################################


## Set up ---------------------------------------------------------------------------------------------------
rm(list=ls())
username <- Sys.getenv("USER")

library(data.table)
library(dplyr)

#System configuration
l_root <- "/FILEPATH/"

#Set directories
read_dir_flywheel <- paste0(l_root, "FILEPATH/")
read_dir <- paste0(l_root, "FILEPATH/")

save_dir <- paste0(l_root, "FILEPATH/")


## Functions ---------------------------------------------------------------------------------------------------
## Create spectral power long data set by time, brain region, waveform for  use across analyses ----------------------------
clean_eeg_power <- function(dt){
  #Subset dataset to spectral power variables and clean (only relevant variables + subject ID + good channels)
  long_power <- copy(as.data.frame(dt))
  
  #Retain power variables  
  keep_strings <- c("Gamma", "gamma", "Theta", "theta", "Alpha", "Beta", "beta", "Delta", "delta", "subject_id", "time", "eeg_age_days")
  long_power = long_power[,grepl(paste(keep_strings, collapse = "|"), x = names(long_power))]
  
  
  #For spectral power, retain average across both hemispheres
  #remove_strings <- c("R_", "L_", "central_", "BL_", "ROI", "_ave")
  remove_strings <- c("ROI")
  long_power = long_power[,!grepl(paste(remove_strings, collapse = "|"), x = names(long_power))]
  
  #Remove unnecessary string
  names(long_power) <- sub("_T1", "", names(long_power))
  names(long_power) <- sub("_ave", "", names(long_power))
  names(long_power) <- sub("BL_", "", names(long_power))
  
  #Average left and right hemisphere power by waveform
  long_power <- long_power %>%
    mutate(across(where(is.character) & !subject_id, as.numeric))
  setDT(long_power)
  long_power[, power_Gamma_front := (power_Gamma_lfront+power_Gamma_rfront)/2]
  long_power[, power_Delta_front := (power_Delta_lfront+power_Delta_rfront)/2]
  long_power[, power_Theta_front := (power_Theta_lfront+power_Theta_rfront)/2]
  long_power[, power_HighAlpha_front := (power_HighAlpha_lfront+power_HighAlpha_rfront)/2]
  long_power[, power_LowAlpha_front := (power_LowAlpha_lfront+power_LowAlpha_rfront)/2]
  long_power[, power_Beta_front := (power_Beta_lfront+power_Beta_rfront)/2]
  
  long_power[, power_Gamma_temp := (power_Gamma_ltemp+power_Gamma_rtemp)/2]
  long_power[, power_Delta_temp := (power_Delta_ltemp+power_Delta_rtemp)/2]
  long_power[, power_Theta_temp := (power_Theta_ltemp+power_Theta_rtemp)/2]
  long_power[, power_HighAlpha_temp := (power_HighAlpha_ltemp+power_HighAlpha_rtemp)/2]
  long_power[, power_LowAlpha_temp := (power_LowAlpha_ltemp+power_LowAlpha_rtemp)/2]
  long_power[, power_Beta_temp := (power_Beta_ltemp+power_Beta_rtemp)/2]
  
  long_power[, power_Gamma_par := (power_Gamma_lpar+power_Gamma_rpar)/2]
  long_power[, power_Delta_par := (power_Delta_lpar+power_Delta_rpar)/2]
  long_power[, power_Theta_par := (power_Theta_lpar+power_Theta_rpar)/2]
  long_power[, power_HighAlpha_par := (power_HighAlpha_lfront+power_HighAlpha_rpar)/2]
  long_power[, power_LowAlpha_par := (power_LowAlpha_lpar+power_LowAlpha_rpar)/2]
  long_power[, power_Beta_par := (power_Beta_lpar+power_Beta_rpar)/2]
  
  long_power <- as.data.frame(long_power)
  remove_strings <- c("rpar", "lpar", "rtemp", "ltemp", "rfront", "lfront")
  long_power = long_power[,!grepl(paste(remove_strings, collapse = "|"), x = names(long_power))]
  
  #Create variables for waveform and brain region
  par_string <- c("subject_id", "eeg_age_days", "time", "_par")
  occip_string <- c("subject_id", "eeg_age_days", "time", "_occip")
  front_string <- c("subject_id", "eeg_age_days", "time", "_front")
  temp_string <- c("subject_id", "eeg_age_days", "time", "_temp")  
  
  par_power = long_power[,grepl(paste(par_string, collapse = "|"), x = names(long_power))]
  occip_power = long_power[,grepl(paste(occip_string, collapse = "|"), x = names(long_power))]
  front_power = long_power[,grepl(paste(front_string, collapse = "|"), x = names(long_power))]
  temp_power = long_power[,grepl(paste(temp_string, collapse = "|"), x = names(long_power))]
  
  par_power$brain_region = "Parietal"
  occip_power$brain_region = "Occipital"
  front_power$brain_region = "Frontal"
  temp_power$brain_region = "Temporal"
  
  names(par_power) <- sub("_par", "", names(par_power))
  names(occip_power) <- sub("_occip", "", names(occip_power))
  names(front_power) <- sub("_front", "", names(front_power))
  names(temp_power) <- sub("_temp", "", names(temp_power))
  
  long_power2 <- plyr::rbind.fill(par_power, occip_power, front_power, temp_power)
  
  high_alpha_string <- c("subject_id", "eeg_age_days", "timepoint", "brain_region", "HighAlpha")
  low_alpha_string <- c("subject_id", "eeg_age_days", "timepoint", "brain_region", "LowAlpha")
  beta_string <- c("subject_id", "eeg_age_days", "timepoint", "brain_region", "Beta")
  gamma_string <- c("subject_id", "eeg_age_days", "timepoint", "brain_region", "Gamma")  
  delta_string <- c("subject_id", "eeg_age_days", "timepoint", "brain_region", "Delta")  
  theta_string <- c("subject_id", "eeg_age_days", "timepoint", "brain_region", "Theta")  
  
  high_alpha_power = long_power2[,grepl(paste(high_alpha_string, collapse = "|"), x = names(long_power2))]
  low_alpha_power = long_power2[,grepl(paste(low_alpha_string, collapse = "|"), x = names(long_power2))]
  beta_power = long_power2[,grepl(paste(beta_string, collapse = "|"), x = names(long_power2))]
  delta_power = long_power2[,grepl(paste(delta_string, collapse = "|"), x = names(long_power2))]
  theta_power = long_power2[,grepl(paste(theta_string, collapse = "|"), x = names(long_power2))]
  gamma_power = long_power2[,grepl(paste(gamma_string, collapse = "|"), x = names(long_power2))]
  
  high_alpha_power$waveform = "High alpha"
  low_alpha_power$waveform = "Low alpha"
  beta_power$waveform = "Beta"
  delta_power$waveform = "Delta"
  gamma_power$waveform = "Gamma"
  theta_power$waveform = "Theta"
  
  
  names(high_alpha_power) <- sub("_HighAlpha", "", names(high_alpha_power))
  names(low_alpha_power) <- sub("_LowAlpha", "", names(low_alpha_power))
  names(beta_power) <- sub("_Beta", "", names(beta_power))
  names(gamma_power) <- sub("_Gamma", "", names(gamma_power))
  names(delta_power) <- sub("_Delta", "", names(delta_power))
  names(theta_power) <- sub("_Theta", "", names(theta_power))
  
  long_power3 <- plyr::rbind.fill(high_alpha_power, low_alpha_power, beta_power, delta_power, theta_power, gamma_power)
  
  #Remove all NA rows for these variables
  long_power4 <- long_power3[!is.na(long_power3$power),] #remove all NA rows
  long_power4 <- long_power4[!is.na(long_power4$eeg_age_days),] #remove all NA rows
  
  long_power4 <- as.data.table(long_power4)
  
  #Update update age variable
  long_power4[, eeg_age_days := round(long_power4$eeg_age_days, digits=0)]
  
  # Identify outliers based on Percent_Good_Chans_Selected being less than 70%
  #long_power4[, is_outlier := 0] 
  #long_power4[, is_outlier := fifelse(Percent_Good_Chans_Selected < 70, 1, 0), by = .(subject_id, time)]
  
  return(long_power4)
}  


## Data manipulation ---------------------------------------------------------------------------------------------------
orig_dt <- as.data.table(read.csv(paste0(read_dir_flywheel, "UNITY-dataTemplate_MRI_clean.csv")))

df_total <- read_excel(paste0(read_dir_flywheel, "SPACE_longitudinal_EEG spectralparameterization.xlsx"))
df_power <- read_excel(paste0(read_dir, "ZAF_GUGULETHU_SPACE_EEG_2024_2025_POWERBAND_LONGITUDINAL_Y2025M08D25.XLSX"))

#Harmonize subject ID and age variables for EEG files
setnames(df_total, "EEG_subject_id", "subject_id")
setnames(df_power, "eeg_age", "eeg_age_days")

df_power <- subset(df_power, select = -c(timepoint))

#More wrangling
##power
df_power <- clean_eeg_power(df_power)
df_power <- rename(df_power, measure = waveform)  
df_power <- rename(df_power, value = power)

theta_beta_wide <- dcast( # calculate theta/beta ratio
  df_power[measure %in% c("Theta", "Beta")],
  subject_id + eeg_age_days + brain_region ~ measure,
  value.var = "value",
  fun.aggregate = mean
)
theta_beta_wide[, value := Theta / Beta]
theta_beta_wide[, measure := "theta_beta_ratio"]
theta_beta_long <- theta_beta_wide[, .(subject_id, eeg_age_days, brain_region, measure, value)]
df_power <- rbind(df_power, theta_beta_long)

alpha_delta_wide <- dcast( # calculate high alpha/delta ratio
  df_power[measure %in% c("High alpha", "Delta")],
  subject_id + eeg_age_days + brain_region ~ measure,
  value.var = "value",
  fun.aggregate = mean
)
alpha_delta_wide[, value := `High alpha` / Delta]
alpha_delta_wide[, measure := "high_alpha_delta_ratio"]
alpha_delta_long <- alpha_delta_wide[, .(subject_id, eeg_age_days, brain_region, measure, value)]
df_power <- rbind(df_power, alpha_delta_long)
df_power <- as.data.table(df_power)

##total power
df_total$Average_aper_exponent <- as.numeric(df_total$Average_aper_exponent)
df_total$Average_aper_offset <- as.numeric(df_total$Average_aper_offset)
df_total <- as.data.table(df_total[,c("subject_id", "eeg_age_days", 
                                      "avg_alpha_by_chan", "Average_aper_exponent", "Average_aper_offset")])
df_total <- df_total[!is.na(eeg_age_days),]
df_total[, avg_alpha_by_chan:= as.numeric(avg_alpha_by_chan)]
df_total <- pivot_longer(df_total, cols = c("avg_alpha_by_chan", "Average_aper_exponent", "Average_aper_offset"), values_to = "value", names_to = "measure") #make long
df_total <- df_total[!is.na(df_total$value),] #remove all NA rows
df_total$brain_region <- "All"

df_total <- as.data.table(df_total)

#Combine
list_combine <- list(df_power, df_total)
df_combined <- rbindlist(list_combine, use.names = TRUE)

## Read out datasets ---------------------------------------------------------------------------------------------------
#Convert to wide
df_combined[brain_region == "All", brain_region := ""]
df_combined[, full_measure := paste0(brain_region, measure)]
df_combined[, full_measure := gsub(" ", "_", full_measure)]

df_combined_fix_dup <- summarise(
  group_by(
    df_combined,
    subject_id,
    eeg_age_days,
    full_measure
  ),
  value = mean(value, na.rm = TRUE)
)

df_combined_wide <- pivot_wider(
  data = df_combined_fix_dup,
  names_from = full_measure,
  values_from = value
)

df_combined_wide <- as.data.table(df_combined_wide) 

#Create binary alpha frequency variable
has_total_subset <- df_combined_wide[!is.na(Average_aper_exponent),]
has_total_subset[, alpha_peak_binary := ifelse(is.na(avg_alpha_by_chan), 0, 1)]

has_total_subset <- subset(has_total_subset, select = c(subject_id, eeg_age_days, alpha_peak_binary))

df_combined_wide <- merge(df_combined_wide, has_total_subset, by = c("subject_id", "eeg_age_days"), all.x = TRUE)

#Fix naming
new_names <- c(
  "FrontalBeta"                      = "power_Beta_front",
  "FrontalDelta"                     = "power_Delta_front",
  "FrontalGamma"                     = "power_Gamma_front",
  "FrontalHigh_alpha"               = "power_HighAlpha_front",
  "FrontalLow_alpha"                = "power_LowAlpha_front",
  "Frontalhigh_alpha_delta_ratio"   = "high_alpha_delta_ratio_front",
  "Frontaltheta_beta_ratio"         = "theta_beta_ratio_front",
  "FrontalTheta"         = "power_Theta_front",
  
  "OccipitalBeta"                   = "power_Beta_occip",
  "OccipitalDelta"                  = "power_Delta_occip",
  "OccipitalGamma"                  = "power_Gamma_occip",
  "OccipitalHigh_alpha"             = "power_HighAlpha_occip",
  "OccipitalLow_alpha"              = "power_LowAlpha_occip",
  "Occipitalhigh_alpha_delta_ratio"= "high_alpha_delta_ratio_occip",
  "Occipitaltheta_beta_ratio"      = "theta_beta_ratio_occip",
  "OccipitalTheta"         = "power_Theta_occip",
  
  "ParietalBeta"                    = "power_Beta_par",
  "ParietalDelta"                   = "power_Delta_par",
  "ParietalGamma"                   = "power_Gamma_par",
  "ParietalHigh_alpha"             = "power_HighAlpha_par",
  "ParietalLow_alpha"              = "power_LowAlpha_par",
  "Parietalhigh_alpha_delta_ratio" = "high_alpha_delta_ratio_par",
  "Parietaltheta_beta_ratio"       = "theta_beta_ratio_par",
  "ParietalTheta"         = "power_Theta_par",
  
  "TemporalBeta"                    = "power_Beta_temporal",
  "TemporalDelta"                   = "power_Delta_temporal",
  "TemporalGamma"                   = "power_Gamma_temporal",
  "TemporalHigh_alpha"             = "power_HighAlpha_temporal",
  "TemporalLow_alpha"              = "power_LowAlpha_temporal",
  "Temporalhigh_alpha_delta_ratio" = "high_alpha_delta_ratio_temporal",
  "Temporaltheta_beta_ratio"       = "theta_beta_ratio_temporal",
  "TemporalTheta"         = "power_Theta_temporal"
)

setnames(df_combined_wide, old = names(new_names), new = new_names)

#Read out full dataset with all rows and all vars
write.csv(df_combined_wide, paste0(save_dir, "space_all_eeg.csv"), row.names = F)

#Subset to selected analysis variables
keep_cols <- c("subject_id", "eeg_age_days", "alpha_peak_binary", "avg_alpha_by_chan",
               "power_Beta_front", "power_Delta_front", "power_Gamma_front", "power_HighAlpha_front", "power_LowAlpha_front", "high_alpha_delta_ratio_front", "theta_beta_ratio_front",
               "power_Beta_occip", "power_Delta_occip", "power_Gamma_occip", "power_HighAlpha_occip", "power_LowAlpha_occip", "high_alpha_delta_ratio_occip", "theta_beta_ratio_occip",
               "power_Beta_par", "power_Delta_par", "power_Gamma_par", "power_HighAlpha_par", "power_LowAlpha_par", "high_alpha_delta_ratio_par", "theta_beta_ratio_par")

df_combined_wide_subset <- df_combined_wide[, ..keep_cols]

#Read out full dataset with all rows and all vars
write.csv(df_combined_wide_subset, paste0(save_dir, "space_subset_all_eeg_selected_vars.csv"), row.names = F)

#Subset to selected analysis variables
length(unique(df_combined_wide_subset$subject_id))
df_combined_wide_subset <- df_combined_wide_subset[df_combined_wide_subset[, .I[which.min(eeg_age_days)], by = subject_id]$V1]
length(unique(df_combined_wide_subset$subject_id))

#Read out subsetted dataset
write.csv(df_combined_wide_subset, paste0(save_dir, "space_subset_earliest_eeg_selected_vars.csv"), row.names = F)







