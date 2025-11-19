# ============================================================================
# XGBoost Modeling Functions for Neurodevelopmental Risk Prediction
# ============================================================================
# This script contains functions for fitting XGBoost models with:
#   - Automatic parameter tuning to prevent overfitting
#   - Missing data imputation using missForest
#   - SHAP value calculation for model interpretability
#   - Cross-validation performance evaluation
#   - Support for imbalanced binary classification tasks
# ============================================================================

# Load required packages
library(xgboost)      # Gradient boosting framework
library(missForest)   # Random forest-based imputation
library(caret)        # Model evaluation and confusion matrices
library(PRROC)        # Precision-Recall curve and AUC-PR calculation
library(fastshap)     # Fast SHAP value computation
library(shapviz)      # SHAP visualization
library(dplyr)        # Data manipulation

# ============================================================================
# MAIN MODELING FUNCTION
# ============================================================================
# Function to fit XGBoost model with missing data imputation
# 
# Parameters:
#   xgb_data: Data frame containing response and predictor variables
#   response: Name of binary outcome variable (factor with 2 levels)
#   covariates: Vector of predictor variable names
#   max_trials: Maximum number of parameter combinations to try (default: 500)
#   seed: Random seed for reproducibility (default: 123)
#
# Returns:
#   List containing fitted model, performance metrics, SHAP values, and more
# ============================================================================
fit_xgboost_miss <- function(xgb_data, response, covariates, max_trials = 500, seed = 123) {
  set.seed(seed)
  resp_levels <- levels(xgb_data[[response]])
  pos_class <- resp_levels[2]
  
  # --- Build modeling frame and impute with missForest ---
  X <- as.data.frame(xgb_data[, covariates, drop = FALSE])
  imp <- missForest::missForest(X, verbose = FALSE)
  X_imp <- imp$ximp
  df_imp <- cbind(xgb_data[[response]], X_imp)
  names(df_imp)[1] <- response
  
  # --- Calculate scale_pos_weight for class balancing ---
  class_counts <- table(df_imp[[response]])
  n_negative <- class_counts[resp_levels[1]]
  n_positive <- class_counts[pos_class]
  scale_pos_weight <- n_negative / n_positive
  
  # --- Prepare data for XGBoost ---
  # Convert response to numeric (0 for first level, 1 for second level)
  y <- as.numeric(df_imp[[response]]) - 1
  
  # Convert covariates to matrix, handling factors
  X <- df_imp[covariates]
  X_matrix <- model.matrix(~ ., data = X)
  X_matrix <- X_matrix[, colnames(X_matrix) != "(Intercept)", drop = FALSE]
  
  # --- Automatic parameter tuning to prevent overfitting ---
  cat("🔍 Automatically tuning XGBoost parameters to minimize overfitting...\n")
  tuning_result <- tune_xgboost_params(xgb_data = df_imp, response = response, 
                                      covariates = covariates, seed = seed, 
                                      max_trials = max_trials)
  
  # Use tuned parameters
  cat("✅ Using tuned parameters (overfitting gap =", round(tuning_result$best_overfitting, 3), ")\n")
  
  nrounds_to_use <- tuning_result$best_params$nrounds
  early_stopping_to_use <- 5
  
  # --- Display final parameters ---
  cat("🎯 Final XGBoost parameters (tuned):\n")
  cat("   max_depth =", tuning_result$best_params$max_depth, "| eta =", tuning_result$best_params$eta,
      "| subsample =", tuning_result$best_params$subsample, "| colsample_bytree =", tuning_result$best_params$colsample_bytree, "\n")
  cat("   min_child_weight =", tuning_result$best_params$min_child_weight, "| alpha =", tuning_result$best_params$alpha,
      "| lambda =", tuning_result$best_params$lambda, "| nrounds =", nrounds_to_use, "\n\n")
  
  # --- Use tuned best full model; compute in-sample metrics here ---
  xgb_fit <- tuning_result$best_full_model
  predictions <- predict(xgb_fit, X_matrix)
  binary_predictions <- factor(ifelse(predictions > 0.5, pos_class, resp_levels[1]), levels = resp_levels)
  metrics_in <- caret::confusionMatrix(binary_predictions,
                                    df_imp[[response]],
                                    positive = pos_class)
  cm_in <- metrics_in$table
  accuracy_in    <- unname(metrics_in$overall["Accuracy"])
  precision_in   <- unname(metrics_in$byClass["Precision"])
  recall_in      <- unname(metrics_in$byClass["Sensitivity"])
  specificity_in <- unname(metrics_in$byClass["Specificity"])
  f1_score_in    <- unname(metrics_in$byClass["F1"])
  auc_pr_in <- tryCatch({
    pr_curve <- pr.curve(scores.class0 = predictions, 
                         weights.class0 = as.numeric(df_imp[[response]] == pos_class),
                         curve = TRUE)
    pr_curve$auc.integral
  }, error = function(e) NA_real_)
  
  # --- CV performance from tuner (pooled across folds) ---
  cv_perf <- tuning_result$cv_performance
  cv_cm <- cv_perf$confusion_matrix_cv
  cv_recall <- cv_perf$recall_cv
  cv_precision <- cv_perf$precision_cv
  cv_f1 <- cv_perf$f1_score_cv
  cv_auc <- cv_perf$auc_pr_cv
  cv_specificity <- cv_perf$specificity_cv
  cv_accuracy <- cv_perf$accuracy_cv
  # SDs no longer used (pooled metrics only)
  
  # --- Variable importance ---
  importance_scores <- xgb.importance(feature_names = colnames(X_matrix), model = xgb_fit)
  var_importance <- data.frame(
    Variable = importance_scores$Feature,
    Gain = importance_scores$Gain,
    Cover = importance_scores$Cover,
    Frequency = importance_scores$Frequency,
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(Gain))
  cat("    Top 10 variables:", paste(head(var_importance$Variable, 10), collapse = ", "), "\n")
  
  # --- SHAP Analysis using fastshap + shapviz ---
  tryCatch({
    # Create explainer object using fastshap
    explainer <- fastshap::explain(
      xgb_fit, 
      X = X_matrix,
      nsim = 100,  # Number of Monte Carlo simulations
      pred_wrapper = function(object, newdata) {
        predict(object, newdata)
      }
    )
    
    # Create shapviz object from fastshap results
    sv <- shapviz(explainer, X = X_matrix)
    
    # Get mean absolute SHAP values for variable importance
    shap_importance <- colMeans(abs(sv$S))
    shap_df <- data.frame(
      Variable = names(shap_importance),
      SHAP_Importance = shap_importance,
      stringsAsFactors = FALSE
    ) %>%
      arrange(desc(SHAP_Importance))
    
    cat("    Top 10 SHAP variables:", paste(head(shap_df$Variable, 10), collapse = ", "), "\n")
    
    # Store shapviz object for plotting
    shap_viz_obj <- sv
  }, error = function(e) {
    cat("    SHAP calculation failed:", e$message, "\n")
    shap_df <- NULL
    shap_viz_obj <- NULL
  })
  
  # --- Console summary ---
  cat("    In-sample Performance | Recall:", round(recall_in, 3),
      "| Precision:", round(precision_in, 3),
      "| F1:", round(f1_score_in, 3),
      "| AUC-PR:", round(as.numeric(auc_pr_in), 3),
      "| Specificity:", round(specificity_in, 3),
      "| Accuracy:", round(accuracy_in, 3), "\n")
  
  cat("    CV Performance (pooled) | Recall:", round(cv_recall, 3),
      "| Precision:", round(cv_precision, 3),
      "| F1:", round(cv_f1, 3),
      "| AUC-PR:", round(cv_auc, 3),
      "| Specificity:", round(cv_specificity, 3),
      "| Accuracy:", round(cv_accuracy, 3), "\n")
  cat("    CV Confusion Matrix (pooled):\n")
  print(cv_cm)
  
  # --- Return ---
  return(list(
    xgb_fit = xgb_fit,
    imputed_data = df_imp,
    X_matrix = X_matrix,
    performance = list(
      # In-sample metrics (ordered by importance for imbalanced data)
      recall_in = recall_in,
      precision_in = precision_in,
      f1_score_in = f1_score_in,
      auc_pr_in = as.numeric(auc_pr_in),
      specificity_in = specificity_in,
      accuracy_in = accuracy_in,
      confusion_matrix_in = cm_in,
      # CV metrics (pooled across folds)
      confusion_matrix_cv = cv_cm,
      recall_cv = cv_recall,
      precision_cv = cv_precision,
      f1_score_cv = cv_f1,
      auc_pr_cv = cv_auc,
      specificity_cv = cv_specificity,
      accuracy_cv = cv_accuracy
    ),
    variable_importance = var_importance,
    shap_importance = shap_df,
    shap_viz_obj = shap_viz_obj,
    predictions = predictions,             # numeric probs for positive class
    binary_predictions = binary_predictions,
    pos_class = pos_class,
    scale_pos_weight = scale_pos_weight,
    tuning_result = tuning_result,        # Results from parameter tuning
    final_params = tuning_result$best_params            # Final tuned parameters
  ))
}

# ============================================================================
# AUTOMATIC PARAMETER TUNING FUNCTION
# ============================================================================
# Function to automatically select optimal XGBoost parameters
# Uses cross-validation to find parameters that minimize overfitting
# while maximizing performance on imbalanced binary classification tasks
tune_xgboost_params <- function(xgb_data, response, covariates, seed = 123, 
                               max_trials = 100) {
  set.seed(seed)
  resp_levels <- levels(xgb_data[[response]])
  pos_class <- resp_levels[2]
  
  # Calculate scale_pos_weight
  class_counts <- table(xgb_data[[response]])
  n_negative <- class_counts[resp_levels[1]]
  n_positive <- class_counts[pos_class]
  scale_pos_weight <- n_negative / n_positive
  
  # Prepare matrices
  y <- as.numeric(xgb_data[[response]]) - 1
  X <- xgb_data[covariates]
  X_matrix <- model.matrix(~ ., data = X)
  X_matrix <- X_matrix[, colnames(X_matrix) != "(Intercept)", drop = FALSE]
  
  # --- Parameter grid for tuning ---
  param_grid <- expand.grid(
    max_depth = c(1, 2, 3),
    eta = c(0.03, 0.05, 0.08, 0.1),
    subsample = c(0.6, 0.7, 0.8),
    colsample_bytree = c(0.6, 0.7, 0.8),
    min_child_weight = c(3, 5, 8),
    alpha = c(0.1, 0.5, 1.0),
    lambda = c(1.0, 2.0, 3.0),
    nrounds = c(50, 100, 150)
  )
  
  # Randomly sample from parameter grid
  set.seed(seed)
  sampled_params <- param_grid[sample(nrow(param_grid), min(max_trials, nrow(param_grid))), ]
  
  best_params <- NULL
  best_overfitting <- Inf
  best_score <- -Inf

  # Holders for artifacts corresponding to current bests
  best_cv_performance <- NULL
  best_full_model <- NULL
  
  for (i in 1:nrow(sampled_params)) {
    params <- sampled_params[i, ]
    
    # Set up parameters
    xgb_params <- list(
      objective = "binary:logistic",
      eval_metric = "aucpr",
      max_depth = params$max_depth,
      eta = params$eta,
      subsample = params$subsample,
      colsample_bytree = params$colsample_bytree,
      min_child_weight = params$min_child_weight,
      scale_pos_weight = scale_pos_weight,
      alpha = params$alpha,
      lambda = params$lambda,
      nthread = 1
    )
    
    # Quick CV evaluation (5-fold for speed)
    cv_folds <- caret::createFolds(xgb_data[[response]], k = 5, list = TRUE, returnTrain = FALSE)
    # Hold pooled predictions and truths across folds
    pooled_pred_prob <- c()
    pooled_pred_class <- character(0)
    pooled_true <- character(0)
    
    for (fold in 1:5) {
      test_indices <- cv_folds[[fold]]
      train_data <- xgb_data[-test_indices, ]
      test_data <- xgb_data[test_indices, ]
      
      y_train <- as.numeric(train_data[[response]]) - 1
      X_train <- model.matrix(~ ., data = train_data[covariates])
      X_train <- X_train[, colnames(X_train) != "(Intercept)", drop = FALSE]
      
      y_test <- as.numeric(test_data[[response]]) - 1
      X_test <- model.matrix(~ ., data = test_data[covariates])
      X_test <- X_test[, colnames(X_test) != "(Intercept)", drop = FALSE]
      
      # Fit model
      xgb_cv <- xgboost(
        data = X_train,
        label = y_train,
        params = xgb_params,
        nrounds = params$nrounds,
        verbose = 0,
        early_stopping_rounds = 5
      )
      
      # Predictions and accumulation for pooling
      test_pred_prob <- predict(xgb_cv, X_test)
      test_pred_class <- ifelse(test_pred_prob > 0.5, pos_class, resp_levels[1])
      pooled_pred_prob <- c(pooled_pred_prob, test_pred_prob)
      pooled_pred_class <- c(pooled_pred_class, test_pred_class)
      pooled_true <- c(pooled_true, as.character(test_data[[response]]))
      
      # No per-fold metrics needed; pooled metrics computed after loop
    }
    
    # Build pooled confusion matrix and metrics
    pooled_pred_factor <- factor(pooled_pred_class, levels = resp_levels)
    pooled_true_factor <- factor(pooled_true, levels = resp_levels)
    pooled_cm <- caret::confusionMatrix(pooled_pred_factor, pooled_true_factor, positive = pos_class)
    pooled_tbl <- pooled_cm$table
    pooled_accuracy <- unname(pooled_cm$overall["Accuracy"])
    pooled_precision <- unname(pooled_cm$byClass["Precision"])
    pooled_recall <- unname(pooled_cm$byClass["Sensitivity"])
    pooled_specificity <- unname(pooled_cm$byClass["Specificity"])
    pooled_f1 <- if (is.na(pooled_recall) || is.na(pooled_precision) || (pooled_recall + pooled_precision) == 0) {
      0
    } else {
      2 * (pooled_recall * pooled_precision) / (pooled_recall + pooled_precision)
    }
    pooled_auc_pr <- tryCatch({
      pr_curve <- pr.curve(scores.class0 = pooled_pred_prob,
                           weights.class0 = as.numeric(pooled_true_factor == pos_class),
                           curve = TRUE)
      pr_curve$auc.integral
    }, error = function(e) NA_real_)
    
    cv_perf_this <- list(
      confusion_matrix_cv = pooled_tbl,
      recall_cv = pooled_recall,
      precision_cv = pooled_precision,
      f1_score_cv = pooled_f1,
      auc_pr_cv = as.numeric(pooled_auc_pr),
      specificity_cv = pooled_specificity,
      accuracy_cv = pooled_accuracy
    )
    
    # Fit full model for in-sample performance
    xgb_full <- xgboost(
      data = X_matrix,
      label = y,
      params = xgb_params,
      nrounds = params$nrounds,
      verbose = 0,
      early_stopping_rounds = 5
    )
    
    # In-sample predictions
    in_sample_pred <- predict(xgb_full, X_matrix)
    auc_pr_in <- tryCatch({
      pr_curve <- pr.curve(scores.class0 = in_sample_pred, 
                           weights.class0 = as.numeric(xgb_data[[response]] == pos_class),
                           curve = TRUE)
      pr_curve$auc.integral
    }, error = function(e) NA_real_)
    # In-sample F1 (using 0.5 threshold)
    in_sample_class <- factor(ifelse(in_sample_pred > 0.5, pos_class, resp_levels[1]), levels = resp_levels)
    in_metrics <- caret::confusionMatrix(in_sample_class, xgb_data[[response]], positive = pos_class)
    f1_in <- unname(in_metrics$byClass["F1"])
    
    # Calculate overfitting gap (F1: in-sample vs pooled CV)
    overfitting_gap <- as.numeric(f1_in) - pooled_f1
    
    # Update best when pooled F1 improves (initialize on first iteration)
    if (is.null(best_params) || pooled_f1 > best_score) {
      best_params <- params
      best_overfitting <- overfitting_gap
      best_score <- pooled_f1
      best_cv_performance <- cv_perf_this
      best_full_model <- xgb_full
    }
  }
  
  return(list(
    best_params = best_params,
    best_score = best_score,
    best_overfitting = best_overfitting,
    best_full_model = best_full_model,
    cv_performance = best_cv_performance
  ))
}

# ============================================================================
# SHAP VISUALIZATION FUNCTION
# ============================================================================
# Function to create SHAP plots
create_shap_page_xgb <- function(rf_results, title = NULL, max_display = NULL) {
  # Define descriptive names for predictors
  predictor_names <- c(
    "sexMale" = "Sex", "age_scan" = "Age at Scan", "age_eeg" = "Age at EEG", "age_cog" = "Age at Score", 
    "meduHigh" = "Maternal Education", "ga" = "Gestation Weeks", 
    "ma" = "Maternal Age", "fs" = "Family Size",
    "bw" = "Birth Weight", "bl" = "Birth Length", 
    "hcmri" = "Head Circumference", 
    "waz" = "WAZ", "haz" = "HAZ", "whz" = "WHZ", "cw" = "Weight", "ch" = "Height",
    "LowAlpha" = "Low Alpha", "HighAlpha" = "High Alpha",
    "Beta" = "Beta", "Gamma" = "Gamma", "ThetaBeta" = "Theta Beta Ratio", 
    "HighAlphaDelta" = "High Alpha Delta Ratio", "AlphaPeak1" = "Alpha Peak",
    "caudate" = "Caudate", "putamen" = "Putamen", "thalamus" = "Thalamus", 
    "ventricles" = "Ventricles", 
    "anterior_callosum" = "Anterior Callosum", 
    "central_callosum" = "Central Callosum",
    "posterior_callosum" = "Posterior Callosum",
    "corpus_callosum" = "Corpus Callosum",
    "supratentorial_csf" = "Supratentorial CSF", 
    "supratentorial_tissue" = "Supratentorial Tissue",
    "cerebellum" = "Cerebellum", "global_pallidus" = "Globus Pallidus", 
    "sc_csf" = "SC CSF", "sc_gray_matter" = "SC Gray Matter", "sc_white_matter" = "SC White Matter", "sc_corpus_callosum" = "SC Corpus Callosum", 
    "sc_caudate" = "SC Caudate", "sc_lentiform" = "SC Lentiform", "sc_hippocampus" = "SC Hippocampus", "sc_thalamus" = "SC Thalamus"#,
    # # Indicator variables for missing values
    # "sex_miss" = "Sex MI", "age_scan_miss" = "Age at Scan MI", "age_eeg_miss" = "Age at EEG MI", 
    # "age_cog_miss" = "Age at Score MI", "medu_miss" = "Maternal Education MI", "ga_miss" = "Gestation Weeks MI",
    # "ma_miss" = "Maternal Age MI", "fs_miss" = "Family Size MI", "bw_miss" = "Birth Weight MI", 
    # "bl_miss" = "Birth Length MI", "hcmri_miss" = "Head Circumference MI", "waz_miss" = "WAZ MI",
    # "haz_miss" = "HAZ MI", "whz_miss" = "WHZ MI", "cw_miss" = "Weight MI", "ch_miss" = "Height MI",
    # "LowAlpha_miss" = "Low Alpha MI", "HighAlpha_miss" = "High Alpha MI", "Beta_miss" = "Beta MI", 
    # "Gamma_miss" = "Gamma MI", "Delta_miss" = "Delta MI", "Ratio_miss" = "High Alpha Delta Ratio MI", 
    # "AlphaPeak_miss" = "Alpha Peak MI", "caudate_miss" = "Caudate MI", "putamen_miss" = "Putamen MI", 
    # "thalamus_miss" = "Thalamus MI", "ventricles_miss" = "Ventricles MI", 
    # "anterior_callosum_miss" = "Anterior Callosum MI", "central_callosum_miss" = "Central Callosum MI",
    # "posterior_callosum_miss" = "Posterior Callosum MI", "supratentorial_csc_miss" = "Supratentorial CSF MI", 
    # "supratentorial_tissue_miss" = "Supratentorial Tissue MI", "cerebellum_miss" = "Cerebellum MI", 
    # "global_pallidus_miss" = "Globus Pallidus MI", "sc_csc_miss" = "SC CSF MI", 
    # "sc_gray_matter_miss" = "SC Gray Matter MI", "sc_white_matter_miss" = "SC White Matter MI", 
    # "sc_caudate_miss" = "SC Caudate MI", "sc_lentiform_miss" = "SC Lentiform MI", 
    # "sc_hippocampus_miss" = "SC Hippocampus MI"
  )
  
  # Helper function to clean and relabel
  relabel <- function(x) {
    x <- gsub("(_residual|_percentile)$", "", x)      # drop suffixes
    ifelse(x %in% names(predictor_names), predictor_names[x], x)
  }
  
  # Helper function to get non-zero SHAP values for a single rf_results
  get_nonzero_shap_vars <- function(rf_result) {
    shap_data <- rf_result$shap_viz_obj$S
    non_zero_vars <- colnames(shap_data)[colSums(abs(shap_data)) > 0]
    return(non_zero_vars)
  }
  
  # Determine if input is single or multiple rf_results
  if (is.list(rf_results) && "shap_viz_obj" %in% names(rf_results)) {
    # Single rf_results case
    if (is.null(max_display)) {
      max_display <- length(get_nonzero_shap_vars(rf_results))
    }
    if (is.null(title)) {
      title <- "GSED"
    }
    
    # Global importance
    global <- sv_importance(rf_results$shap_viz_obj, max_display = max_display) + 
      ggtitle(paste0(title, " - SHAP")) +
      scale_y_discrete(labels = relabel)
    
    return(list(global))
    
  } else if (is.list(rf_results) && length(rf_results) == 4) {
    # Four rf_results case
    if (is.null(title)) {
      title <- c("BSID CCS", "BSID LCS", "BSID MCS", "GSED")
    }
    
    # Get union of non-zero SHAP variables across all four models
    if (is.null(max_display)) {
      all_nonzero_vars <- unique(unlist(lapply(rf_results, get_nonzero_shap_vars)))
      max_display <- length(all_nonzero_vars)
    }
    
    # Create plots for each model
    plots <- list()
    for (i in 1:4) {
      plots[[i]] <- sv_importance(rf_results[[i]]$shap_viz_obj, max_display = max_display) + 
        ggtitle(paste0(title[i], " - SHAP")) +
        scale_y_discrete(labels = relabel)
    }
    
    return(plots)
  } else {
    stop("rf_results must be either a single rf_results object or a list of 4 rf_results objects")
  }
  
  # # Local beeswarm
  # beeswarm <- sv_importance(rf_results$shap_viz_obj, kind = "beeswarm", max_display = max_display) + 
  #   ggtitle(paste0(title, " - Local SHAP")) +
  #   scale_y_discrete(labels = relabel)
  # 
  # for (i in seq_along(beeswarm$layers)) {
  #   if (inherits(beeswarm$layers[[i]]$geom, "GeomPoint")) {
  #     beeswarm$layers[[i]]$aes_params$alpha <- 0.8
  #     beeswarm$layers[[i]]$aes_params$size <- 0.5
  #   }
  # }
  # 
  # list(global, beeswarm)
}

# ============================================================================
# PERFORMANCE SUMMARY HELPER FUNCTION
# ============================================================================
# Helper function to convert model results into a one-row performance summary
# Used for creating performance comparison tables
# ============================================================================
perf_row <- function(name, res) {
  p <- res$performance
  tibble::tibble(
    Model = name,
    # `Recall (in)`      = round(p$recall_in, 3),
    # `Precision (in)`   = round(p$precision_in, 3),
    # `F1 (in)`          = round(p$f1_score_in, 3),
    # `AUC-PR (in)`      = round(p$auc_pr_in, 3),
    # `Specificity (in)` = round(p$specificity_in, 3),
    # `Accuracy (in)`    = round(p$accuracy_in, 3),
    # 
    # `Recall (CV)`      = sprintf("%.3f ± %.3f", p$cv_recall,      p$cv_recall_sd),
    # `Precision (CV)`   = sprintf("%.3f ± %.3f", p$cv_precision,   p$cv_precision_sd),
    # `F1 (CV)`          = sprintf("%.3f ± %.3f", p$cv_f1_score,    p$cv_f1_score_sd),
    # `AUC-PR (CV)`      = sprintf("%.3f ± %.3f", p$cv_auc_pr,      p$cv_auc_pr_sd),
    # `Specificity (CV)` = sprintf("%.3f ± %.3f", p$cv_specificity, p$cv_specificity_sd),
    # `Accuracy (CV)`    = sprintf("%.3f ± %.3f", p$cv_accuracy,    p$cv_accuracy_sd)
    `Recall (CV)`      = round(p$recall_cv, 3),
    `Precision (CV)`   = round(p$precision_cv, 3),
    `F1 (CV)`          = round(p$f1_score_cv, 3),
    `AUC-PR (CV)`      = round(p$auc_pr_cv, 3),
    `Specificity (CV)` = round(p$specificity_cv, 3),
    `Accuracy (CV)`    = round(p$accuracy_cv, 3)
  )
}
