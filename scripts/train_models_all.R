# =============================================================
# train_models_all.R
# =============================================================

#' Trains and evaluates classification models on global methylation matrices.
#'
#' This function uses full methylation matrices (all probes per panel) to train
#' supervised models for disease classification.
#' Supports SVM and Random Forest (ranger). Optionally includes Age/Sex as covariates.
#'
#' @param root_dir Path to the directory containing global train/test CSVs
#' @param prob_threshold Probability cutoff to assign predicted classes (default = 0.5)
#' @param output_file CSV file to save model evaluation results
#' @param type ML algorithm: "svmLinear", "svmRadial" or "ranger"
#' @param valid_split_vars Character vector of extra covariates to include
#' @param panel_name Name of the gene panel used (for logging)
#' @param seed Random seed for reproducibility
#' @param tune_length Number of hyperparameter combinations to test in caret::train
#' @param contrast_suffix Suffix used to identify clinical contrasts
#' @param sample_train_df Optional sample sheet for training set (required if using covariates)
#' @param sample_test_df Optional sample sheet for test set (required if using covariates)
#'
#' @return A list with two elements: results (data.frame with metrics), and models (list of trained caret objects)
train_models_all <- function(
    root_dir = "global",
    prob_threshold = 0.5,
    output_file = "ml_results_all.csv",
    type = "svmLinear",
    valid_split_vars = character(),
    panel_name = "default_panel",
    seed = 123,
    tune_length = 5,
    contrast_suffix = "vsCTRL",
    sample_train_df = NULL,
    sample_test_df  = NULL
) {
  set.seed(seed)
  cat("Training GLOBAL ML models for panel:", panel_name, "...\n")
  
  library(caret)
  library(dplyr)
  
  results_df <- data.frame()
  model_list <- list()
  
  # Ensure Sample_Name is used as rownames in sample sheets
  if (!is.null(sample_train_df) && "Sample_Name" %in% colnames(sample_train_df)) {
    rownames(sample_train_df) <- sample_train_df$Sample_Name
  }
  if (!is.null(sample_test_df) && "Sample_Name" %in% colnames(sample_test_df)) {
    rownames(sample_test_df) <- sample_test_df$Sample_Name
  }
  
  # Define cross-validation strategy
  tr_control <- caret::trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE,
    savePredictions = "final",
    summaryFunction = twoClassSummary
  )
  
  # Recursively locate global training files
  all_train_files <- list.files(root_dir, pattern = paste0("^train_all_.*", contrast_suffix, "\\.csv$"), recursive = TRUE, full.names = TRUE)
  
  for (train_file in all_train_files) {
    contrast <- sub("^train_all_", "", basename(train_file))
    contrast <- sub("\\.csv$", "", contrast)
    
    test_file <- file.path(dirname(train_file), paste0("test_all_", contrast, ".csv"))
    
    message("Processing global contrast: ", contrast)
    
    # Read train matrix
    train_df <- tryCatch({
      df <- read.csv(train_file, row.names = 1)
      if (ncol(df) < 2 || nrow(df) == 0) stop("Invalid train file.")
      df
    }, error = function(e) {
      message("[ERROR] Reading train file: ", train_file)
      return(NULL)
    })
    
    # Read test matrix
    test_df <- tryCatch({
      df <- read.csv(test_file, row.names = 1)
      if (ncol(df) < 2 || nrow(df) == 0) stop("Invalid test file.")
      df
    }, error = function(e) {
      message("[ERROR] Reading test file: ", test_file)
      return(NULL)
    })
    
    if (is.null(train_df) || is.null(test_df)) next
    
    # Determine positive label
    positive_label <- sub("vsCTRL$", "", contrast)
    
    if (positive_label == "neuro") {
      train_df$Sample_Group <- ifelse(grepl("^CTRL", rownames(train_df)), "CTRL", "neuro")
      test_df$Sample_Group  <- ifelse(grepl("^CTRL", rownames(test_df)),  "CTRL", "neuro")
    } else {
      keep_regex <- paste0("^(CTRL|", positive_label, ")")
      train_df <- train_df[grepl(keep_regex, rownames(train_df)), , drop = FALSE]
      test_df  <- test_df [grepl(keep_regex, rownames(test_df )), , drop = FALSE]
      train_df$Sample_Group <- ifelse(grepl(paste0("^", positive_label), rownames(train_df)), positive_label, "CTRL")
      test_df$Sample_Group  <- ifelse(grepl(paste0("^", positive_label), rownames(test_df)),  positive_label, "CTRL")
    }
    
    # Set factor levels
    train_df$Sample_Group <- factor(train_df$Sample_Group, levels = c("CTRL", positive_label))
    test_df$Sample_Group  <- factor(test_df$Sample_Group,  levels = c("CTRL", positive_label))
    
    # Add Age/Sex from sample sheets if ranger is used
    if (type == "ranger" && length(valid_split_vars) > 0) {
      if (is.null(sample_train_df) || is.null(sample_test_df)) {
        stop("You must provide sample_train_df and sample_test_df to use valid_split_vars.")
      }
      
      extra_train <- sample_train_df %>%
        dplyr::filter(Sample_Name %in% rownames(train_df)) %>%
        dplyr::select(all_of(c("Sample_Name", valid_split_vars))) %>%
        tibble::remove_rownames() %>%
        tibble::column_to_rownames("Sample_Name") %>%
        .[rownames(train_df), , drop = FALSE]
      
      extra_test <- sample_test_df %>%
        dplyr::filter(Sample_Name %in% rownames(test_df)) %>%
        dplyr::select(all_of(c("Sample_Name", valid_split_vars))) %>%
        tibble::remove_rownames() %>%
        tibble::column_to_rownames("Sample_Name") %>%
        .[rownames(test_df), , drop = FALSE]
      
      for (var in valid_split_vars) {
        if (var == "Sex") {
          extra_train[[var]] <- as.factor(extra_train[[var]])
          extra_test[[var]]  <- as.factor(extra_test[[var]])
        } else if (var == "Age") {
          extra_train[[var]] <- as.numeric(extra_train[[var]])
          extra_test[[var]]  <- as.numeric(extra_test[[var]])
        }
      }
      
      train_df <- cbind(train_df, extra_train)
      test_df  <- cbind(test_df,  extra_test)
    }
    
    # Skip if only one class is present
    if (length(unique(train_df$Sample_Group)) < 2) {
      message("[SKIP] Only one class in training set for ", contrast)
      next
    }
    
    # Train ML model
    model_fit <- tryCatch({
      if (type == "ranger") {
        predictors <- setdiff(colnames(train_df), "Sample_Group")
        formula_obj <- as.formula(paste("Sample_Group ~", paste(predictors, collapse = " + ")))
        train(formula_obj,
              data = train_df,
              method = "ranger",
              trControl = tr_control,
              tuneLength = tune_length)
      } else {
        train(Sample_Group ~ .,
              data = train_df,
              method = type,
              trControl = tr_control,
              tuneLength = tune_length)
      }
    }, error = function(e) {
      message("[ERROR] Training failed in ", contrast)
      return(NULL)
    })
    
    if (is.null(model_fit)) next
    
    # Predict on test set
    pred_prob  <- predict(model_fit, newdata = test_df, type = "prob")
    pred_class <- ifelse(pred_prob[[positive_label]] > prob_threshold, positive_label, "CTRL")
    pred_class <- factor(pred_class, levels = c("CTRL", positive_label))
    
    # Evaluate predictions
    cm <- confusionMatrix(pred_class, test_df$Sample_Group, positive = positive_label)
    
    row_res <- data.frame(
      comparison        = contrast,
      region            = "GLOBAL",
      p_value           = unname(cm$overall["AccuracyPValue"]),
      mcnemar_p_value   = unname(cm$overall["McnemarPValue"]),
      NIR               = unname(cm$overall["AccuracyNull"]),
      accuracy          = unname(cm$overall["Accuracy"]),
      balanced_accuracy = unname(cm$byClass["Balanced Accuracy"]),
      kappa             = unname(cm$overall["Kappa"]),
      stringsAsFactors  = FALSE
    )
    
    # Store results
    results_df <- bind_rows(results_df, row_res)
    model_list[[contrast]] <- model_fit
  }
  
  # Warn if no models were trained
  if (nrow(results_df) == 0) {
    warning("No models were trained. Check input files or structure.")
    return(NULL)
  }
  
  # Save results to file
  write.csv(results_df, output_file, row.names = FALSE)
  cat("ML GLOBAL results saved to:", output_file, "\n\n")
  
  invisible(list(results = results_df, models = model_list))
}