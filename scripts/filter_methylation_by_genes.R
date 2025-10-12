# =============================================================
# filter_methylation_by_genes.R
# =============================================================

#' Filters EPIC methylation probes based on a user-defined gene list.
#'
#' This function retrieves gene coordinates from Ensembl (via biomaRt),
#' formats them as BED, and intersects them with EPIC probe coordinates
#' using bedtools. It returns a filtered methylation matrix and summary stats.
#'
#' @param gene_list Vector of gene symbols
#' @param probes_annotation_file CSV file with EPIC probe coordinates
#' @param methylation_file CSV methylation matrix (rows = probes, cols = samples)
#' @param output_file Path to save the filtered methylation matrix
#' @param tmp_genes_file Path to save temporary BED file for genes
#' @param tmp_probes_file Path to save temporary BED file for probes
#' @param tmp_intersect_file Path to save bedtools intersection output
#' @param panel_name Name of the panel (used for logging)
#'
#' @return Summary stats
filter_methylation_by_genes <- function(
    gene_list,
    probes_annotation_file,
    methylation_file,
    output_file,
    tmp_genes_file,
    tmp_probes_file,
    tmp_intersect_file,
    panel_name
) {
  cat("\n Filtering probes based on provided gene list...\n")
  
  # Count initial number of genes
  n_genes_extracted <- length(gene_list)
  
  # Get genomic coordinates for each gene from Ensembl
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  genes_coords <- getBM(
    attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype"),
    filters = "hgnc_symbol",
    values = gene_list,
    mart = ensembl
  ) %>%
    dplyr::filter(chromosome_name %in% 1:22) %>%
    dplyr::mutate(chromosome_name = paste0("chr", chromosome_name))
  
  # Count mapped genes
  n_genes_mapped <- genes_coords %>%
    dplyr::pull(hgnc_symbol) %>%
    unique() %>%
    length()
  
  # Format gene coordinates as BED
  genes_bed <- genes_coords %>%
    dplyr::select(
      chrom = chromosome_name,
      start = start_position,
      end = end_position,
      name = hgnc_symbol
    ) %>%
    dplyr::mutate(
      start = as.integer(start),
      end   = as.integer(end)
    ) %>%
    dplyr::arrange(chrom, start)
  
  # Format EPIC probe annotations as BED
  probes <- read.csv(probes_annotation_file)
  epic_bed <- probes %>%
    dplyr::select(chrom = chr, start = pos, name = Name) %>%
    dplyr::mutate(
      chrom  = ifelse(grepl("^chr", chrom), chrom, paste0("chr", chrom)),
      start0 = as.integer(start) - 1L,
      end    = as.integer(start)
    ) %>%
    dplyr::select(chrom, start = start0, end, name) %>%
    dplyr::arrange(chrom, start)
  
  # Write temporary BED files for bedtools
  write_tsv(genes_bed, tmp_genes_file, col_names = FALSE)
  write_tsv(epic_bed,  tmp_probes_file, col_names = FALSE)
  
  # Intersect using bedtools
  bedtools_cmd <- sprintf(
    "bedtools intersect -a %s -b %s -wa -wb > %s",
    tmp_probes_file, tmp_genes_file, tmp_intersect_file
  )
  system(bedtools_cmd)
  
  if (!file.exists(tmp_intersect_file)) stop("Error: Intersection not generated")
  
  # Extract intersected probes
  intersection <- readr::read_tsv(
    tmp_intersect_file,
    col_names = c("chr", "start", "end", "probe_id",
                  "gene_chr", "gene_start", "gene_end", "gene_name"),
    show_col_types = FALSE
  )
  intersected_probes <- unique(intersection$probe_id)
  
  # Filter methylation matrix by probe IDs
  methylation <- read.csv(methylation_file, row.names = 1)
  filtered_methylation <- methylation[rownames(methylation) %in% intersected_probes, ]
  
  # Save filtered methylation matrix
  write.csv(filtered_methylation, output_file, row.names = TRUE)
  
  # Log basic stats
  cat("Filtering completed.\n")
  cat("Extracted genes:", n_genes_extracted, "\t Mapped genes:", n_genes_mapped, "\n")
  cat("Original probes:", nrow(methylation), "\t Filtered probes:", nrow(filtered_methylation), "\n\n")
  
  # Return stats for tracking
  return(list(
    total_genes_input    = n_genes_extracted,
    genes_mapped         = n_genes_mapped,
    probes_total_filtered = nrow(filtered_methylation)
  ))
}