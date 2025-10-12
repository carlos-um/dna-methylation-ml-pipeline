# =============================================================
# train_models_loci.R
# =============================================================

#' Trains and evaluates classification models for each genomic locus.
#'
#' For each locus and clinical contrast, this function loads
#' locus-specific methylation matrices and trains supervised classifiers using caret.
#' Supports SVM and Random Forest (ranger). Optionally includes Age/Sex as covariates.
#'
#' @param root_dir Directory containing locus subfolders (each with contrast folders)
#' @param prob_threshold Probability cutoff for assigning binary classes (default = 0.5)
#' @param output_file CSV file where evaluation metrics will be saved
#' @param type ML algorithm to use: "svmLinear", "svmRadial", or "ranger"
#' @param valid_split_vars Character vector of additional covariates
#' @param panel_name Name of the gene panel (used for logging)
#' @param seed Random seed
#' @param tune_length Number of hyperparameter combinations to try
#' @param contrast_suffix Pattern used to identify contrast folders
#' @param sample_train_df Optional sample metadata for training (used if covariates included)
#' @param sample_test_df Optional sample metadata for testing (used if covariates included)
#'
#' @return A list with two elements: results (data.frame with metrics), and models (list of trained caret objects)
train_models_loci <- function(
    root_dir = "loci",                
    prob_threshold = 0.5,          
    output_file = "ml_results.csv",  
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
  cat("Training ML models by loci for", panel_name, "...\n")
  
  library(caret)
  library(dplyr)
  
  results_df <- data.frame()  # Data frame to store evaluation results
  model_list <- list()        # List to store trained models
  
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
  
  # List all functional loci
  regions <- list.dirs(root_dir, full.names = FALSE, recursive = FALSE)
  
  for (region in regions) {
    comparison_dirs <- list.dirs(file.path(root_dir, region), full.names = FALSE, recursive = FALSE)
    
    for (comp in comparison_dirs) {
      # Skip if the folder name doesn't match the expected contrast suffix
      if (!grepl(paste0(contrast_suffix, "$"), comp)) next
      
      message("Processing: ", region, ";", comp)
      
      # Build file paths to train/test data
      train_file <- file.path(root_dir, region, comp, sprintf("%s_train_ml_%s.csv", region, comp))
      test_file  <- file.path(root_dir, region, comp, sprintf("%s_test_ml_%s.csv",  region, comp))
      
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
      positive_label <- sub("vsCTRL$", "", comp)
      negative_label <- "CTRL"
      levels_vec <- c(negative_label, negative_label)
      
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
        message("[SKIP] Only one class in training set for ", region, "|", comp)
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
        message("[ERROR] Training failed in ", region, "|", comp)
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
        comparison        = comp,
        region            = region,
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
      model_list[[paste(region, comp, sep = "_")]] <- model_fit
    }
  }
  
  # Exit if no models were trained
  if (nrow(results_df) == 0) {
    warning("No models were trained. Check input directories and file structure.")
    return(NULL)
  }
  
  # Save results to file
  write.csv(results_df, output_file, row.names = FALSE)
  cat("ML results saved to:", output_file, "\n\n")
  
  invisible(list(results = results_df, models = model_list))
}
