# =============================================================
# extract_gene_panels.R
# =============================================================

# ---- Function 1: Extract from full TSV export (PanelApp) ----

#' Extracts "Expert Review Green" genes from a Genomics England TSV file.
#' These genes are considered clinically relevant based on expert consensus.
#'
#' @param genomics_england_file Path to the TSV file exported from PanelApp
#' @return A character vector with unique gene symbols
get_green_genes_from_genomics_england <- function(genomics_england_file) {
  cat("Extracting 'Expert Review Green' genes from:", genomics_england_file, "\n")
  
  df_genes <- readr::read_tsv(genomics_england_file, show_col_types = FALSE)
  
  genes_of_interest <- df_genes %>%
    dplyr::filter(grepl("Expert Review Green", `Sources(; separated)`)) %>%
    dplyr::distinct(`Gene Symbol`, .keep_all = TRUE) %>%
    dplyr::pull(`Gene Symbol`) %>%
    na.omit() %>%
    unique()
  
  cat("Extracted", length(genes_of_interest), "unique green genes.\n\n")
  return(genes_of_interest)
}


# ---- Function 2: Extract from simplified CSV file ----

#' Extracts genes marked as "HighEvidence" from a simplified CSV panel file.
#' The file must contain 'GeneSymbol' and 'LevelOfConfidence' columns.
#'
#' @param csv_file Path to the simplified CSV file
#' @return A character vector with unique gene symbols
get_high_evidence_genes_from_csv <- function(csv_file) {
  cat("Extracting 'HighEvidence' genes from:", csv_file, "\n")
  
  df <- readr::read_csv(csv_file, show_col_types = FALSE)
  
  genes_of_interest <- df %>%
    dplyr::filter(LevelOfConfidence == "HighEvidence") %>%
    dplyr::distinct(GeneSymbol, .keep_all = TRUE) %>%
    dplyr::pull(GeneSymbol) %>%
    na.omit() %>%
    unique()
  
  cat("Extracted", length(genes_of_interest), "unique high-evidence genes.\n\n")
  return(genes_of_interest)
}
