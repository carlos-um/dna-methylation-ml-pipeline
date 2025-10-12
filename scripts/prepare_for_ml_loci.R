# =============================================================
# prepare_for_ml_loci.R
# =============================================================

#' Prepares train/test matrices for supervised learning, segmented by genomic loci.
#' For each functional region (promoters, UTRs...), this function creates matched
#' training and test matrices and assigns binary labels per clinical contrast.
#'
#' @param sample_sheet_train_file Metadata for training samples
#' @param sample_sheet_test_file Metadata for test samples
#' @param test_file Global test methylation matrix (probes x samples)
#' @param loci_dir Directory containing locus-specific train matrices
#' @param contrasts Vector of contrasts
#' @param loci Functional loci to process 
#' @return No return value. Saves train/test matrices per locus and contrast.
prepare_for_ml_loci <- function(
    sample_sheet_train_file,
    sample_sheet_test_file,
    test_file,
    loci_dir,
    contrasts = c("DLBvsCTRL", "PDvsCTRL", "PDDvsCTRL", "neurovsCTRL"),
    loci = c("promoter", "tss200", "body", "stExon", "utr5", "tss1500", "utr3", "exonbnd")
) {
  cat("Preparing train/test matrices for ML (segmented by locus)...\n")
  
  # Read sample sheets
  sample_train <- read.csv(sample_sheet_train_file)
  sample_test  <- read.csv(sample_sheet_test_file)
  rownames(sample_train) <- sample_train$Sample_Name
  rownames(sample_test)  <- sample_test$Sample_Name
  
  # Read and prepare global test matrix
  test_mat <- read.csv(test_file, row.names = 1)
  test_df <- as.data.frame(t(test_mat))  # samples as rows
  colnames(test_df) <- rownames(test_mat)  # ensure probe names as column names
  
  # basic validation
  stopifnot(all(rownames(test_df) %in% rownames(sample_test)))
  
  # Iterate over each functional locus
  for (locus in loci) {
    train_file <- file.path(loci_dir, locus, paste0(locus, "_train.csv"))
    if (!file.exists(train_file)) {
      warning("Missing training file for locus: ", locus)
      next
    }
    
    # Read locus-specific training matrix
    train_df <- read.csv(train_file, row.names = 1)
    
    # Validate that all samples exist in the sample sheet
    missing_train <- setdiff(rownames(train_df), rownames(sample_train))
    if (length(missing_train) > 0) {
      warning("Missing sample(s) in sample_train: ", paste(missing_train, collapse = ", "))
      train_df <- train_df[rownames(train_df) %in% rownames(sample_train), ]
    }
    
    # Subset test matrix to probes available in train
    probes_to_use <- colnames(train_df)
    test_df_common <- test_df[, probes_to_use, drop = FALSE]
    
    # Iterate over each clinical contrast
    for (contrast in contrasts) {
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
      test_ids  <- intersect(rownames(test_df_common), rownames(sample_test)[sample_test$Group %in% c(pos, neg)])
      
      X_train <- train_df[train_ids, , drop = FALSE]
      X_test  <- test_df_common[test_ids, , drop = FALSE]
      
      X_train <- X_train[order(rownames(X_train)), , drop = FALSE]
      X_test  <- X_test[order(rownames(X_test)), , drop = FALSE]
      
      sample_train_sub <- sample_train[rownames(X_train), , drop = FALSE]
      sample_test_sub  <- sample_test[rownames(X_test),  , drop = FALSE]
      
      # Add class labels
      X_train$Sample_Group <- ifelse(sample_train_sub$Group %in% pos, pos[1], "CTRL")
      X_test$Sample_Group  <- ifelse(sample_test_sub$Group  %in% pos, pos[1], "CTRL")
      
      X_train$Sample_Group <- factor(X_train$Sample_Group, levels = c("CTRL", pos[1]))
      X_test$Sample_Group  <- factor(X_test$Sample_Group,  levels = c("CTRL", pos[1]))
      
      # Basic validation before saving
      stopifnot(!is.null(rownames(X_test)))
      stopifnot(!is.null(colnames(X_test)))
      
      # Export to files
      contrast_dir <- file.path(loci_dir, locus, contrast)
      if (!dir.exists(contrast_dir)) dir.create(contrast_dir, recursive = TRUE)
      
      write.csv(X_train, file.path(contrast_dir, paste0(locus, "_train_ml_", contrast, ".csv")), row.names = TRUE)
      write.csv(X_test,  file.path(contrast_dir, paste0(locus, "_test_ml_",  contrast, ".csv")),  row.names = TRUE)
    }
  }
  
  cat("All train/test matrices successfully prepared.\n\n")
}
