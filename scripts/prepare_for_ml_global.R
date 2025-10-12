# =============================================================
# prepare_for_ml_global.R
# =============================================================

#' Prepares global train/test matrices for supervised machine learning.
#' This function ensures consistent features (CpG probes) between datasets
#' and assigns binary labels for specified clinical contrasts.
#'
#' @param filtered_train_file Path to the filtered training methylation matrix (probes x samples)
#' @param test_file Path to the global test methylation matrix (probes x samples)
#' @param sample_sheet_train_file Metadata file for training samples (Sample_Name, Group)
#' @param sample_sheet_test_file Metadata file for test samples (Sample_Name, Group)
#' @param output_dir Output directory for saving train/test CSV files
#' @param contrasts Character vector of clinical contrasts
#' @return No return value. Writes train/test CSV files to disk.
prepare_for_ml_global <- function(
    filtered_train_file,
    test_file,
    sample_sheet_train_file,
    sample_sheet_test_file,
    output_dir,
    contrasts = c("DLBvsCTRL", "PDvsCTRL", "PDDvsCTRL", "neurovsCTRL")
) {
  cat("Preparing global train/test matrices for ML (train_all_ / test_all_)...\n")
  
  # Reading methylation matrices
  # Loaded as probes x samples
  train_mat <- read.csv(filtered_train_file, row.names = 1)
  test_mat  <- read.csv(test_file, row.names = 1)

  # Transpose to match ML input format
  train_df <- as.data.frame(t(train_mat))
  test_df  <- as.data.frame(t(test_mat))
  
  # Reading metadata tables (sample sheets)
  sample_train <- read.csv(sample_sheet_train_file)
  sample_test  <- read.csv(sample_sheet_test_file)
  
  # Assign sample names as rownames for later matching
  rownames(sample_train) <- sample_train$Sample_Name
  rownames(sample_test)  <- sample_test$Sample_Name
  
  # Ensure train and test matrices share the same probes
  probes_to_use <- intersect(colnames(train_df), colnames(test_df))
  train_df <- train_df[, probes_to_use, drop = FALSE]
  test_df  <- test_df[, probes_to_use, drop = FALSE]
  
  # Iterate over each clinical (binary) contrast
  for (contrast in contrasts) {
    
    # Define positive and negative groups (disease vs control)
    if (contrast == "DLBvsCTRL") {
      pos <- "DLB"; neg <- "CTRL"
    } else if (contrast == "PDvsCTRL") {
      pos <- "PD"; neg <- "CTRL"
    } else if (contrast == "PDDvsCTRL") {
      pos <- "PDD"; neg <- "CTRL"
    } else if (contrast == "neurovsCTRL") {
      pos <- c("DLB", "PD", "PDD"); neg <- "CTRL"
    } else {
      stop("Unknown contrast: ", contrast)
    }
    
    # Filter and match samples
    train_ids <- intersect(rownames(train_df), rownames(sample_train)[sample_train$Group %in% c(pos, neg)])
    test_ids  <- intersect(rownames(test_df),  rownames(sample_test)[sample_test$Group  %in% c(pos, neg)])
    
    X_train <- train_df[train_ids, , drop = FALSE]
    X_test  <- test_df[test_ids, , drop = FALSE]
    
    X_train <- X_train[order(rownames(X_train)), , drop = FALSE]
    X_test  <- X_test[order(rownames(X_test)), , drop = FALSE]
    
    sample_train_sub <- sample_train[rownames(X_train), , drop = FALSE]
    sample_test_sub  <- sample_test[rownames(X_test),  , drop = FALSE]
    
    # Add class labels
    X_train$Sample_Group <- ifelse(sample_train_sub$Group %in% pos, pos[1], "CTRL")
    X_test$Sample_Group  <- ifelse(sample_test_sub$Group  %in% pos, pos[1], "CTRL")
    
    # Define factor levels with order CTRL < POS
    X_train$Sample_Group <- factor(X_train$Sample_Group, levels = c("CTRL", pos[1]))
    X_test$Sample_Group  <- factor(X_test$Sample_Group,  levels = c("CTRL", pos[1]))
    
    # Basic validation before saving
    stopifnot(!is.null(rownames(X_test)))
    stopifnot(!is.null(colnames(X_test)))
    
    # Export to files
    contrast_dir <- file.path(output_dir, contrast)
    if (!dir.exists(contrast_dir)) dir.create(contrast_dir, recursive = TRUE)
    
    write.csv(X_train, file.path(contrast_dir, paste0("train_all_", contrast, ".csv")), row.names = TRUE)
    write.csv(X_test,  file.path(contrast_dir, paste0("test_all_",  contrast, ".csv")),  row.names = TRUE)
  }
  
  cat("Global train/test matrices generated with prefix 'train_all_' and 'test_all_'.\n\n")
}
