# =============================================================
# split_methylation_by_loci.R
# =============================================================

#' Splits a filtered DNA methylation matrix by functional loci.
#'
#' This function categorizes CpG probes from the Illumina EPIC array into 
#' genomic regions based on UCSC and regulatory annotations. It creates 
#' separate methylation matrices for each region (promoter, TSS200, 
#' gene body, UTRs...) and exports them as CSV files for downstream analysis.
#' Additionally, it saves probe ID lists for Venn diagrams and tracks 
#' summary statistics such as the number of probes and genes per region.
#'
#' @param methylation_file Path to the filtered methylation matrix (CSV)
#' @param output_dir Output directory where region-specific CSVs will be saved
#' @param panel_name Name of the gene panel used (used in file naming)
#' @param venn_output_dir Directory to store probe lists for Venn diagrams
#' @param seed Random seed for reproducibility
#' @param total_genes_input Number of input genes before filtering (for tracking)
#' @param genes_mapped Number of genes successfully mapped to coordinates
#' @param probes_total_filtered Total number of probes after gene-based filtering
#' @return No return value. Outputs CSVs and summary statistics to disk.
split_methylation_by_loci <- function(
    methylation_file = "filtered_methylation.csv",
    output_dir = "loci",
    panel_name = "default_panel", 
    venn_output_dir = "venn_inputs",
    seed = 123,
    total_genes_input = NA,
    genes_mapped = NA,
    probes_total_filtered = NA
) {
  cat("Splitting probes by functional loci...\n")
  set.seed(seed)
  
  # Load methylation matrix and transpose
  methylation_df <- read.csv(methylation_file, row.names = 1)
  methylation_df <- as.data.frame(t(methylation_df))
  
  # Extract probe annotations from EPIC array
  anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  probe_ids <- colnames(methylation_df)
  probe_info <- anno[rownames(anno) %in% probe_ids,
                     c("UCSC_RefGene_Group", "Regulatory_Feature_Group", "UCSC_RefGene_Name")]
  
  # Define promoter probes (TSS + regulatory annotations)
  is_in_tss <- function(group) {
    terms <- unlist(strsplit(group, ";"))
    any(terms %in% c("TSS200", "TSS1500"))
  }
  
  tss_probes <- rownames(probe_info)[sapply(probe_info$UCSC_RefGene_Group, is_in_tss)]
  
  regulatory_criteria <- c("Promoter_Associated", "Promoter_Associated_Cell_type_specific")
  promoter_associated_probes <- rownames(probe_info)[
    grep(paste(regulatory_criteria, collapse = "|"), probe_info$Regulatory_Feature_Group)
  ]
  
  promoter_probes_names <- unique(c(tss_probes, promoter_associated_probes))
  promoter_probes <- probe_info[rownames(probe_info) %in% promoter_probes_names, ]
  
  # Create subsets for each UCSC_RefGene_Group value
  unique_values <- unique(unlist(strsplit(as.character(probe_info$UCSC_RefGene_Group), ";")))
  lista_subsets <- lapply(unique_values, function(val) {
    indices <- sapply(strsplit(as.character(probe_info$UCSC_RefGene_Group), ";"), function(x) val %in% x)
    probe_info[indices, ]
  })
  names(lista_subsets) <- unique_values
  
  # Define output directories by region
  loci_dirs <- list(
    promoter = file.path(output_dir, "promoter"),
    tss200   = file.path(output_dir, "tss200"),
    body     = file.path(output_dir, "body"),
    stExon   = file.path(output_dir, "stExon"),
    utr5     = file.path(output_dir, "utr5"),
    tss1500  = file.path(output_dir, "tss1500"),
    utr3     = file.path(output_dir, "utr3"),
    exonbnd  = file.path(output_dir, "exonbnd")
  )
  for (dir in loci_dirs) if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  
  # Export methylation matrix for each region
  write_subset <- function(subset_name, probe_names) {
    subset_df <- methylation_df[, colnames(methylation_df) %in% probe_names, drop = FALSE]
    write.csv(subset_df,
              file = file.path(loci_dirs[[subset_name]], paste0(subset_name, "_train.csv")),
              row.names = TRUE)
  }
  
  if (!is.null(promoter_probes))           write_subset("promoter", rownames(promoter_probes))
  if (!is.null(lista_subsets$TSS200))      write_subset("tss200",   rownames(lista_subsets$TSS200))
  if (!is.null(lista_subsets$Body))        write_subset("body",     rownames(lista_subsets$Body))
  if (!is.null(lista_subsets$`1stExon`))   write_subset("stExon",   rownames(lista_subsets$`1stExon`))
  if (!is.null(lista_subsets$`5'UTR`))     write_subset("utr5",     rownames(lista_subsets$`5'UTR`))
  if (!is.null(lista_subsets$TSS1500))     write_subset("tss1500",  rownames(lista_subsets$TSS1500))
  if (!is.null(lista_subsets$`3'UTR`))     write_subset("utr3",     rownames(lista_subsets$`3'UTR`))
  if (!is.null(lista_subsets$ExonBnd))     write_subset("exonbnd",  rownames(lista_subsets$ExonBnd))
  
  # Export probe IDs by region for Venn diagram
  for (region_name in names(loci_dirs)) {
    probes_file <- file.path(loci_dirs[[region_name]], paste0(region_name, "_train.csv"))
    if (file.exists(probes_file)) {
      df <- read.csv(probes_file, row.names = 1)
      probe_ids <- colnames(df)
      
      region_dir <- file.path(venn_output_dir, region_name)
      if (!dir.exists(region_dir)) dir.create(region_dir, recursive = TRUE)
      
      out_path <- file.path(region_dir, paste0(panel_name, "_probes.csv"))
      write.csv(data.frame(probe_id = probe_ids), out_path, row.names = FALSE)
    }
  }
  
  # Compute region-wise statistics (n_probes, n_genes)
  loci_stats <- list()
  for (region_name in names(loci_dirs)) {
    probes_file <- file.path(loci_dirs[[region_name]], paste0(region_name, "_train.csv"))
    if (file.exists(probes_file)) {
      df <- read.csv(probes_file, row.names = 1)
      probe_ids <- colnames(df)
      genes <- unique(unlist(strsplit(as.character(probe_info[probe_ids, "UCSC_RefGene_Name"]), ";")))
      loci_stats[[region_name]] <- list(
        n_probes = length(probe_ids),
        n_genes  = length(genes)
      )
    }
  }
  
  # Save region statistics summary file
  save_summary_statistics(
    panel_name            = panel_name,
    output_dir            = dirname(output_dir),
    total_genes_input     = total_genes_input,
    genes_mapped          = genes_mapped,
    probes_total_filtered = probes_total_filtered,
    loci_stats            = loci_stats
  )
  
  cat("Locus-based split completed for:", panel_name, "\n\n")
}