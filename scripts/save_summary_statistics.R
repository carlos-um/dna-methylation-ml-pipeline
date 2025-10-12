# =============================================================
# save_summary_statistics.R
# =============================================================

#' Saves statistics for a given gene panel.
#'
#' This function compiles and exports a summary table that includes the number
#' of input genes, the number of genes successfully mapped to genomic 
#' coordinates, and the number of methylation probes retained after filtering.
#' When available, it also includes per-locus probe and gene counts. The 
#' summary is saved as a CSV file for later inspection and reproducibility.
#'
#' @param panel_name Name of the gene panel being processed
#' @param output_dir Output directory to store the summary CSV
#' @param total_genes_input Total number of gene symbols input to filtering
#' @param genes_mapped Number of genes successfully mapped to coordinates
#' @param probes_total_filtered Total number of probes retained
#' @param loci_stats Named list with per-region stats (n_probes, n_genes)
#' @return No return value. Writes a summary CSV file to disk.
save_summary_statistics <- function(
    panel_name,
    output_dir,
    total_genes_input,
    genes_mapped,
    probes_total_filtered,
    loci_stats = NULL
) {
  # Define path for summary CSV file
  summary_file <- file.path(output_dir, paste0("summary_statistics_", panel_name, ".csv"))
  
  # Base summary: global statistics
  summary_df <- data.frame(
    panel_name            = panel_name,
    total_genes_input     = total_genes_input,
    genes_mapped          = genes_mapped,
    probes_total_filtered = probes_total_filtered,
    stringsAsFactors = FALSE
  )
  
  # Optional: add per-locus statistics
  if (!is.null(loci_stats)) {
    for (region in names(loci_stats)) {
      summary_df[[paste0(region, "_n_probes")]] <- loci_stats[[region]]$n_probes
      summary_df[[paste0(region, "_n_genes")]]  <- loci_stats[[region]]$n_genes
    }
  }
  
  # Export summary table as CSV
  write.csv(summary_df, summary_file, row.names = FALSE)
  cat("Summary saved to:", summary_file, "\n")
}