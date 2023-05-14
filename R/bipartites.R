#' Compute links between TFs and DNA regions
#'
#' Return list of links between peaks and TFs, based on their binding motifs
#' locations on a reference genome.
#' Currently based on Signac AddMotifs function (--> motifmachR, itself based on
#' MOODs algorithm).
#'
#' @param tfs vector(character) - List of tfs considered.
#' @param peaks vector(character) - List of peaks.
#' @param peak_sep1 (character) - Separator between chromosme number
#' and starting coordinates (e.g. : chr1/100_1000 --> sep1='/').
#' @param peak_sep2 (character) - Separator between chromosme number
#' and starting coordinates (e.g. : chr1/100_1000 --> sep1='_').
#' @param genome (BSGenome) - Genome sequences on which motifs positions
#' will be searched for.
#' @param gene.range (gene.range object) - TO DO.
#' @param motifs  (PWMatrixList) List of PWMatrix (PWMs of motifs).
#' @param tf2motifs (data.frame) Corresponding table between PWMatrix names
#' and binding TFs.
#' @param store_bipartite (bool) - Save the bipartite directly
#' (\code{TRUE}, default) or return without saving on disk (\code{FALSE}).
#' @param output_file (character) - Name of the output_file
#' (if store_bipartite == \code{TRUE}).
#' @param verbose (integer) Display function messages.
#' Set to 0 for no message displayed, >= 1 for more details.
#'
#' @return (data.frame) Return list of the links betweeen TFs and peaks.
#' @export
#'
#' @examples TO DO. Same than UNIT test.
bipartite_tfs2peaks <- function(
  hummus_object,
  tf_expr_assay = "RNA",
  peak_assay = "peaks",
  tf_network_name = NULL,
  peak_network_name = NULL,
  genome,
  store_bipartite = FALSE,
  output_file = NULL,
  verbose = 1) {

  if (verbose > 0) {
    cat("Computing TF-peak bipartite\n")
  }
  # Check if tf_gene_assay is NULL
  if (!is.null(tf_expr_assay)) {
    # Check if the gene assay is present in the seurat object
    if (!tf_expr_assay %in% names(hummus_object@assays)) {
      stop("The gene assay is not present in the seurat object")
    }
    # Get TFs expressed in  assay AND having known binding motifs
    tfs_use <- get_tfs(hummus_object,
                       assay = tf_expr_assay,
                       store_tfs = FALSE,
                       verbose=verbose)
  }
  else { # No filtering on expression assay, use all TFs
    tfs_use <- unique(hummus_object@motifs_db@tf2motifs$tf)
  }

  # Check if the peak assay is present in the seurat object
  if (!peak_assay %in% names(hummus_object@assays)) {
    stop("The peak assay is not present in the seurat object")
  }
  # Check if the peak assay is a ChromatinAssay object
  if (!inherits(hummus_object@assays[[peak_assay]],
                     "ChromatinAssay")) {
    stop("The peak assay is not a ChromatinAssay object 
    or does not have annotations (gene.range object))")
  }
  # Check if the peak assay has gene.range annotations
  if (is.null(Signac::Annotation(hummus_object[[peak_assay]]))) {
      stop("The peak assay does not have annotations (gene.range object)")
  }


  motif_pos <- Signac::AddMotifs(
    object = hummus_object[[peak_assay]],
    genome = genome,
    pfm = hummus_object@motifs_db@motifs #add verbose options
  )

  ## The 17 following lines are inspired from the Pando package :
  # https://github.com/quadbiolab/Pando/blob/main/R/regions.R
  # Add TF info for motifs
  if (verbose > 0) {
    cat("\tAdding TF info\n")
  }

  # Spread dataframe to sparse matrix
  tf2motifs <- hummus_object@motifs_db@tf2motifs
  tf2motifs <- dplyr::'%>%'(tf2motifs, dplyr::select("motif" = 1, "tf" = 2)) # Select motif and tf columns
  tf2motifs <- dplyr::'%>%'(tf2motifs, dplyr::distinct()) # Remove duplicates
  tf2motifs <- dplyr::'%>%'(tf2motifs, dplyr::mutate(val = 1)) # Add value column
  tf2motifs <- dplyr::'%>%'(tf2motifs, # Spread TFs
                  tidyr::pivot_wider(names_from = "tf",
                                     values_from = val,
                                     values_fill = 0)
                            )
  tf2motifs <- dplyr::'%>%'(tf2motifs, tibble::column_to_rownames("motif")) # Set motif as rownames
  tf2motifs <- dplyr::'%>%'(tf2motifs, as.matrix()) # Convert to matrix
  tf2motifs <- dplyr::'%>%'(tf2motifs, Matrix::Matrix(sparse = TRUE)) # Convert to sparse matrix

  if (length(tfs_use) == 0) { # If no TFs are found in the dataset
    stop("None of the provided TFs were found in the dataset.
    Consider providing a custom motif-to-TF map as `motif_tfs`")
  }

  # Get TF peak links
  TFs_Peaks <- motif_pos@motifs@data %*% tf2motifs[, tfs_use]
  # TFs_Peaks <- TFs_Peaks[, colnames(TFs_Peaks) %in% tfs]

   # Keep only the TFs that are in our scRNA-seq dataset
  tfs2peaks <- expand.grid(rownames(TFs_Peaks),
                           colnames(TFs_Peaks))[as.vector(TFs_Peaks > 0), ]
                          # TF-peak links
  colnames(tfs2peaks) <- c("peak", "TF")     # set column names

  # Save TF-peak links
  if (store_bipartite) {
    if (is.null(output_file)) {
      stop("Please provide an output file name")
    }
    write.table(tfs2peaks,
                output_file,
                col.names = TRUE,
                row.names = FALSE,
                quote = FALSE,
                sep = "\t")
  }
  if (verbose > 0) {
    cat("\tReturning TF-peak links as bipartite object\n")
  }
  # Return TF-peak links

  # Set default names for the networks if not provided
  if (is.null(tf_network_name)) {
    tf_network_name <- tf_expr_assay
  }
  if (is.null(peak_network_name)) {
    peak_network_name <- peak_assay
  }

  bipartite_tf_peak <- new("bipartite",
                           "network" = tfs2peaks,
                           "multiplex_left" = peak_network_name,
                           "multiplex_right" = tf_network_name)
  return(bipartite_tf_peak) # Return TF-peak bipartite object
}




#' Compute links between DNA regions and genenames
#'
#' Return list of links between peaks and genes,
#' based on the distance between peaks and gene's TSS location
#' from gene.range annotations.
#' Call find_peaks_near_genes function, that can use different methods.
#'
#' @param genes vector(character) - List of genes.
#' @param peaks vector(character) - List of peaks.
#' @param peak_sep1 (character) - Separator between chromosme number
#'  and starting coordinates (e.g. : chr1/100_1000 --> sep1='/').
#' @param peak_sep2 (character) - Separator between chromosme number
#'  and starting coordinates (e.g. : chr1/100_1000 --> sep1='_').
#' @param gene.range (gene.range object) - Gene range object.
#' @param peak_to_gene_method (character) - Method to map peaks to near gene
#' * \code{'Signac'} - Signac method (default).
#' * \code{'GREAT'} - not implemented yet.
#' @param upstream (int) - size of the window upstream the TSS considered
#' @param downstream (int) - size of the window downstream the TSS considered
#' @param only_tss (bool) - Associated peaks in the window size from TSS only (\code{TRUE}) or aroud the whole gene body (\code{FALSE}).
#' @param store_bipartite (bool) - Save the bipartite directly (\code{TRUE}, default) or return without saving on disk (\code{FALSE}).
#' @param output_file (character) - Name of the output_file (if store_bipartite == \code{TRUE}).
#'
#' @return (data.frame) Return list of the links betweeen peaks and genes.
#' @export
#'
#' @examples
bipartite_peaks2genes <- function(
  seurat_object,
  gene_assay = 'RNA',
  peak_assay = 'peaks',
  gene_network_name = NULL,
  peak_network_name = NULL,
  peak_to_gene_method = "Signac",
  upstream = 500,
  downstream = 500,
  only_tss = TRUE,
  store_bipartite = FALSE,
  output_file = NULL) {
  # Check if the gene assay is present in the seurat object
  if (!gene_assay %in% names(seurat_object@assays)) {
    stop("The gene assay is not present in the seurat object")
  }
  # Check if the peak assay is present in the seurat object
  else if (!peak_assay %in% names(seurat_object@assays)) {
    stop("The peak assay is not present in the seurat object")
  }
  # Check if the peak assay is a ChromatinAssay object
  else if (!inherits(seurat_object@assays[[peak_assay]],
                     "ChromatinAssay")) {
    stop("The peak assay is not a ChromatinAssay object 
    or does not have annotations (gene.range object))")
  }
  # Check if the peak assay has gene.range annotations
  else if (is.null(Signac::Annotation(seurat_object[[peak_assay]]))) {
      stop("The peak assay does not have annotations (gene.range object)")
  }

  # Find candidate regions near gene bodies
  peaks_near_genes <- find_peaks_near_genes(
                        peaks = seurat_object[[peak_assay]]@ranges,
                        method = peak_to_gene_method,
                        genes = Signac::Annotation(seurat_object[[peak_assay]]),
                        upstream = upstream,
                        downstream = downstream,
                        only_tss = only_tss)
  # Aggregate candidate regions to gene bodies (peak to gene matrix)
  peaks2genes <- aggregate_matrix(Matrix::t(peaks_near_genes),
                                 groups = colnames(peaks_near_genes),
                                 fun = "sum")
  # Keep only the genes that are in our scRNA-seq dataset
  peaks2genes <- peaks2genes[rownames(peaks2genes) 
                 %in% rownames(seurat_object@assays[[gene_assay]]), ]
  # Remove rows/cols with only zeros
  peaks2genes <- peaks2genes[Matrix::rowSums(peaks2genes) != 0,
                             Matrix::colSums(peaks2genes) != 0]
  # peak-gene links
  peaks2genes <- expand.grid(rownames(peaks2genes),
                             colnames(peaks2genes))[as.vector(peaks2genes > 0), ]
  colnames(peaks2genes) <- c("gene", "peak") # set column names

  if (store_bipartite) {
    if (is.null(output_file)) {
      stop("Please provide an output file name")
    }
    write.table(peaks2genes,
                output_file,
                col.names = TRUE,
                row.names = FALSE,
                quote = FALSE,
                sep = "\t")
  }
  # Set default names for the networks if not provided
  if (is.null(gene_network_name)) {
    gene_network_name <- gene_assay
  }
  if (is.null(peak_network_name)) {
    peak_network_name <- peak_assay
  }

  # Return atac-rna bipartite
  bipartite_atac_rna <- new("bipartite",
                           "network" = peaks2genes,
                           "multiplex_left" = gene_network_name,
                           "multiplex_right" = peak_network_name)
  return(bipartite_atac_rna)
}

#' @title Associate peaks to genes based on distance to TSS (or gene body)
#'
#' @param peaks vector(character) - List of peaks.
#' @param genes vector(character) - List of genes.
#' @param sep vector(character) - Separator between chromosome,
#'            start and end position. Default: c('-', '-').
#' @param method (character) - Method to use. Default: "Signac".
#' * \code{'Signac'} - Use Signac::Extend to extend genes.
#' * \code{'GREAT'} - Not implemented yet.
#' @param upstream (int) - Upstream distance from TSS
#' to consider as potential promoter.
#' @param downstream (int) - Downstream distance from TSS
#' to consider as potential promoter.
#' @param extend (int) - Integer defining the distance from the upstream
#' and downstream of the basal regulatory region. Used only by method 'GREAT'.
#' @param only_tss (logical) - If TRUE, only TSS will be considered.
#' @param verbose (logical) - If TRUE, print progress messages.
#'
#' @return (matrix) - Matrix of peaks x genes with 1 if peak is near gene.
#' @export
#'
#' @examples TODO
find_peaks_near_genes <- function(
  peaks,
  genes,
  sep = c("-", "-"),
  method = c("Signac", "GREAT"),
  upstream = 100000,
  downstream = 0,
  extend = 1000000,
  only_tss = FALSE,
  verbose = TRUE
) {
  # Match arg
  method <- match.arg(method)

  if (method == "Signac") {

    if (only_tss) {
      genes <- IRanges::resize(x = genes, width = 1, fix = "start")
    }
    genes_extended <- suppressWarnings(
      expr = Signac::Extend(
        genes, upstream = upstream, downstream = downstream
      )
    )
    overlaps <- IRanges::findOverlaps(
      query = peaks,
      subject = genes_extended,
      type = "any",
      select = "all"
    )
    hit_matrix <- Matrix::sparseMatrix(
      i = S4Vectors::queryHits(overlaps),
      j = S4Vectors::subjectHits(overlaps),
      x = 1,
      dims = c(length(peaks), length(genes_extended))
    )
    rownames(hit_matrix) <- Signac::GRangesToString(grange = peaks, sep = sep)
    colnames(hit_matrix) <- genes_extended$gene_name

  } else if (method == "_____GREAT") {

    # Read gene annotation (Ensembl v93, GRCh38)
    utils::data(EnsDb.Hsapiens.v93.annot.UCSC.hg38, envir = environment())
    gene_annot_use <- EnsDb.Hsapiens.v93.annot.UCSC.hg38[
      which(EnsDb.Hsapiens.v93.annot.UCSC.hg38$gene_name %in% genes$gene_name),
    ]
    gene_annot_tss <- select(as_tibble(gene_annot_use),
                             seqnames, "start" = tss, "end" = tss, strand)

    # Create GRanges object storing the TSS information
    tss <- GRanges(gene_annot_use)

    # Define basal regulatory region (promoter region)
    # as 5 kb upstream + 1 kb downstream of the TSS
    basal_reg <- suppressWarnings(
      expr = Signac::Extend(
        tss, upstream = upstream, downstream = downstream
      )
    )

    # Step 1 - get peaks overlap with basal regulatory region
    basal_overlaps <- suppressWarnings(IRanges::findOverlaps(
      query = peaks,
      subject = basal_reg,
      type = "any",
      select = "all",
      minoverlap = 2
    ))

    peak_all <- Signac::GRangesToString(grange = peaks, sep = sep)
    basal_peak_mapped_idx <- queryHits(basal_overlaps)
    basal_mapped_peaks <- unique(peak_all[basal_peak_mapped_idx])
    n1 <- length(basal_mapped_peaks)

    # Step 2: for the peaks not overlapped with basal regulatory regions,
    # check whether they located within gene body of any genes
    peak_unmapped_idx <- setdiff(seq(length(peak_all)), basal_peak_mapped_idx)
    peak_unmapped <- peak_all[peak_unmapped_idx]
    peak_unmapped_region <- Signac::StringToGRanges(peak_unmapped)

    # Create GRanges object storing annotated gene boundary
    gene_bound <- GRanges(gene_annot_use)
    body_overlaps <- IRanges::findOverlaps(
      query = peak_unmapped_region,
      subject = gene_bound,
      type = "any",
      select = "all",
      minoverlap = 2
    )
    body_peak_mapped_idx <- peak_unmapped_idx[queryHits(body_overlaps)]
    body_mapped_peaks <- unique(peak_all[body_peak_mapped_idx])
    n2 <- length(body_mapped_peaks)
    peak_mapped_idx <- c(basal_peak_mapped_idx, body_peak_mapped_idx)

    # Step 3: for the peaks not overlapped with regulatory regions of any genes,
    # check whether they overlap with extended regulatory region. 
    # i.e. +/- 1MB of basal regulatory region
    peak_unmapped_idx <- setdiff(seq(length(peak_all)), peak_mapped_idx)
    peak_unmapped <- peak_all[peak_unmapped_idx]
    peak_unmapped_region <- Signac::StringToGRanges(peak_unmapped)
    extend_reg <- suppressWarnings(
      expr = Signac::Extend(
        basal_reg, upstream = extend, downstream = extend
      )
    )

    # Get overlap between unmapped_peak_region and extended regulatory region
    extended_overlaps <- suppressWarnings(IRanges::findOverlaps(
      query = peak_unmapped_region,
      subject = extend_reg,
      type = "any",
      select = "all",
      minoverlap = 2
    ))
    extended_peak_mapped_idx <- peak_unmapped_idx[queryHits(extended_overlaps)]
    extended_mapped_peaks <- unique(peak_all[extended_peak_mapped_idx])
    n3 <- length(extended_mapped_peaks)

    hit_matrix <- Matrix::sparseMatrix(
      i = c(basal_peak_mapped_idx,
            body_peak_mapped_idx,
            extended_peak_mapped_idx),
      j = c(subjectHits(basal_overlaps),
            subjectHits(body_overlaps),
            subjectHits(extended_overlaps)),
      x = 1,
      dims = c(length(peaks), length(basal_reg))
    )
    rownames(hit_matrix) <- peak_all
    colnames(hit_matrix) <- c(basal_reg$gene_name)
  } else {
    stop("method must be either 'Signac' or 'GREAT' ; 
          please check that current version of HuMMuS
          already accepts GREAT as a method.")
  }
  return(hit_matrix)
}


#' @title Filter peaks to those overlapping specific (regulatory) elements
#' @description Function to reduce list of "Peaks" to the ones overlapping with
#' list of "RegEl", e.g. regulatory elements, evolutionary conserved regions
#' 
#' @param Peaks (character) vector of genomic coordinates of peaks
#' @param RegEl (character) vector of genomic coordinates of regulatory elements
#' @param sep_Peak1 (character) separator between chromosome and
#'                              start position of peak
#' @param sep_Peak2 (character) separator between start position
#'                              and end position of peak
#' @param sep_RegEl1 (character) separator between chromosome and
#'                               start position of regulatory element
#' @param sep_RegEl2 (character) separator between start position and
#'                               end position of regulatory element
#'
#' @return (character) vector of genomic coordinates of peaks overlapping
#' @export
#'
#' @examples peaks_in_regulatory_elements(peaks, RegEl)
peaks_in_regulatory_elements <- function(
  Peaks, 
  RegEl, 
  sep_Peak1="-",
  sep_Peak2="-",
  sep_RegEl1="-",
  sep_RegEl2="-"
) {
  # Make sure Peaks and RegEl are unique
  Peaks <- unique(Peaks)
  RegEl <- unique(RegEl)

  # convert genomic corrdinate string to GRanges object
  Peak_GRangesObj <- StringToGRanges(Peaks, sep = c(sep_Peak1, sep_Peak2))
  RegEl_GRangesObj <- StringToGRanges(RegEl, sep = c(sep_RegEl1, sep_RegEl2))

  # find overlap between peaks and regulatory elements
  PeakOverlaps <- IRanges::findOverlaps(query = RegEl_GRangesObj,
                                        subject = Peak_GRangesObj)

  # return peaks that overlapped with regulatory element
  return(Peaks[unique(as.matrix(PeakOverlaps)[, 2])])
}