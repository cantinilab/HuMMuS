#' Compute links between TFs and DNA regions (ATAC peaks)
#'
#' Compute and add bipartite between TFs and DNA regions to hummus object.
#' Links are computed based on the binding motifs of TFs and their locations
#' on a reference genome.
#' Currently based on Signac AddMotifs function (--> motifmachR, itself based on
#' MOODs algorithm).
#'
#' @param hummus_object (hummus_object) - Hummus object.
#' @param tf_expr_assay (character) - Name of assay containing the TF expression
#' data. If NULL, all TFs with a motif are used. Default: "RNA".
#' @param peak_assay (character) - Name of the assay containing the DNA regions
#' (ATAC peaks). Default: "peaks".
#' @param tf_multiplex_name (character) - Name of multiplex containing the TFs.
#' If NULL, the name of the TF assay is used.
#' @param peak_multiplex_name (character) - Name of the multiplex containing the
#' DNA regions (ATAC peaks). If NULL, the name of the peak assay is used.
#' @param genome (BSgenome object) - Reference genome.
#' @param store_network (bool) - Save the bipartite directly
#' (\code{TRUE}, default) or return without saving on disk (\code{FALSE}).
#' @param output_file (character) - Name of the output_file
#' (if store_bipartite == \code{TRUE}). Default: NULL.
#' @param verbose (integer) Display function messages.
#' Set to 0 for no message displayed, >= 1 for more details. Default: 1.
#' @param bipartite_name (character) - Name of bipartite. Default: "tf_peak".
#'
#' @return hummus_object (hummus_object) - Hummus object with TF-peak bipartite
#' added to the multilayer slot
#' @export
#'
#' @examples hummus <- bipartite_tfs2peaks(
#'                      hummus_object = hummus,
#'                      tf_expr_assay = "RNA",
#'                      peak_assay = "peaks",
#'                      tf_multiplex_name = "TF",
#'                      peak_multiplex_name = "peaks",
#'           genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
#'                      store_network = FALSE,
#'                      verbose = 1,
#'                      bipartite_name = "tf_peak")

bipartite_tfs2peaks <- function(
  hummus_object,
  tf_expr_assay = "RNA",
  peak_assay = "peaks",
  tf_multiplex_name = NULL,
  peak_multiplex_name = NULL,
  genome,
  store_network = FALSE,
  output_file = NULL,
  verbose = 1,
  bipartite_name = "tf_peak"
  ) {

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
                       verbose = verbose)
  } else { # No filtering on expression assay, use all TFs with a motif
    if (verbose > 0) {
      cat("No filtering on expression assay, using all TFs with a motif.\n")
    }
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

  # Add motifs to the peaks
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
  # Select motif and tf columns
  tf2motifs <- dplyr::"%>%"(tf2motifs, dplyr::select("motif" = 1, "tf" = 2))
  tf2motifs <- dplyr::"%>%"(tf2motifs, dplyr::distinct()) # Remove duplicates
  # Add value column
  tf2motifs <- dplyr::"%>%"(tf2motifs, dplyr::mutate(val = 1))
  tf2motifs <- dplyr::"%>%"(tf2motifs, # Spread TFs
                  tidyr::pivot_wider(names_from = "tf",
                                     values_from = val,
                                     values_fill = 0)
                            )
  # Set motif as rownames
  tf2motifs <- dplyr::"%>%"(tf2motifs, tibble::column_to_rownames("motif"))
  tf2motifs <- dplyr::"%>%"(tf2motifs, as.matrix()) # Convert to matrix

  # Convert to sparse matrix
  tf2motifs <- dplyr::"%>%"(tf2motifs, Matrix::Matrix(sparse = TRUE))

  if (length(tfs_use) == 0) { # If no TFs are found in the dataset
    stop("None of the provided TFs were found in the dataset.
    Consider providing a custom motif-to-TF map as `motif_tfs`")
  }

  # Get TF peak links
  # Keep only the TFs that are in our tf list
  TFs_Peaks <- motif_pos@motifs@data %*% tf2motifs[, tfs_use]

  # Remove values equal to 0
  tfs2peaks <- expand.grid(rownames(TFs_Peaks),
                           colnames(TFs_Peaks))[as.vector(TFs_Peaks > 0), ]
                          # TF-peak links
  colnames(tfs2peaks) <- c("peak", "TF")     # set column names

  # Save TF-peak links
  store_network(network = tfs2peaks,
                store_network = store_network,
                output_file = output_file,
                verbose = verbose)

  if (verbose > 0) {
    cat("\tReturning TF-peak links as bipartite object\n")
  }

  # Set default names for the networks if not provided
  if (is.null(tf_multiplex_name)) {
    cat("no TF layer name provided, using tf_expr_assay name\n")
    tf_multiplex_name <- tf_expr_assay
  }
  if (is.null(peak_multiplex_name)) {
    cat("no peak layer name provided, using peak_assay name\n")
    peak_multiplex_name <- peak_assay
  }

  # Return tf-peak bipartite
  hummus_object@multilayer@bipartites[[bipartite_name]] <- new("bipartite",
                           "network" = tfs2peaks,
                           "multiplex_left" = peak_multiplex_name,
                           "multiplex_right" = tf_multiplex_name)
  return(hummus_object) # Return TF-peak bipartite object
}




#' Compute links between DNA regions and genenames
#'
#' Compute and add bipartite between DNA regions and genenames to hummus object.
#' Links are computed based on the distance between peaks and gene's TSS
#' location from gene.range annotations.
#' Call find_peaks_near_genes function, that can use different methods.
#'
#' @param hummus_object (hummus_object) - Hummus object.
#' @param gene_assay (character) - Name of assay containing the gene expression
#' data. Default: "RNA".
#' @param peak_assay (character) - Name of the assay containing the DNA regions
#' (ATAC peaks). Default: "peaks".
#' @param gene_multiplex_name (character) - Name of the multiplex containing the
#' genes.
#' If NULL, the name of the gene assay is used.
#' @param peak_multiplex_name (character) - Name of the multiplex containing the
#' DNA regions (ATAC peaks). If NULL, the name of the peak assay is used.
#' @param peak_to_gene_method (character) - Method to use to compute the links
#' between peaks and genes. Default: "Signac".
#' * \code{'Signac'} - Use Signac::Extend to extend genes.
#' * \code{'GREAT'} - Not implemented yet.
#' @param upstream (int) - Upstream distance from TSS
#' to consider as potential promoter.
#' @param downstream (int) - Downstream distance from TSS
#' to consider as potential promoter.
#' @param only_tss (logical) - If TRUE, only TSS will be considered.
#' @param store_network (bool) - Save the bipartite directly
#' (\code{TRUE}, default) or return without saving on disk (\code{FALSE}).
#' @param output_file (character) - Name of the output_file
#' (if store_bipartite == \code{TRUE}). Default: NULL.
#' @param verbose (integer) Display function messages.
#' Set to 0 for no message displayed, >= 1 for more details. Default: 1.
#' @param bipartite_name (character) - Name of bipartite. Default: "atac_rna".
#' 
#' @return hummus_object (hummus_object) - Hummus object w/ atac-rna bipartite
#' added to the multilayer slot
#' @export
#'
#' @examples hummus <- bipartite_peaks2genes(
#'                        hummus_object = hummus,
#'                        gene_assay = "RNA",
#'                        peak_assay = "peaks",
#'                        gene_multiplex_name = "RNA",
#'                        peak_multiplex_name = "peaks",
#'                        peak_to_gene_method = "Signac",
#'                        upstream = 500,
#'                        downstream = 500,
#'                        only_tss = TRUE,
#'                        store_network = FALSE,
#'                        bipartite_name = "atac_rna")

bipartite_peaks2genes <- function(
  hummus_object,
  gene_assay = "RNA",
  peak_assay = "peaks",
  gene_multiplex_name = NULL,
  peak_multiplex_name = NULL,
  peak_to_gene_method = "Signac",
  upstream = 500,
  downstream = 500,
  only_tss = TRUE,
  store_network = FALSE,
  output_file = NULL,
  bipartite_name = "atac_rna"
  ) {
  # Check if the gene assay is present in the hummus object
  if (!gene_assay %in% names(hummus_object@assays)) {
      stop("The gene assay is not present in the hummus object")
  } else if (!peak_assay %in% names(hummus_object@assays)) {
      # Check if the peak assay is present in the hummus object
      stop("The peak assay is not present in the hummus object")
  } else if (!inherits(hummus_object@assays[[peak_assay]],
                      "ChromatinAssay")) {
      # Check if the peak assay is a ChromatinAssay object
      stop("The peak assay is not a ChromatinAssay object 
      or does not have annotations (gene.range object))")
  } else if (is.null(Signac::Annotation(hummus_object[[peak_assay]]))) {
      # Check if the peak assay has gene.range annotations
      stop("The peak assay does not have annotations (gene.range object)")
  }

  # Find candidate regions near gene bodies
  peaks_near_genes <- find_peaks_near_genes(
                        peaks = hummus_object[[peak_assay]]@ranges,
                        method = peak_to_gene_method,
                        genes = Signac::Annotation(hummus_object[[peak_assay]]),
                        upstream = upstream,
                        downstream = downstream,
                        only_tss = only_tss)
  # Aggregate candidate regions to gene bodies (peak to gene matrix)
  peaks2genes <- aggregate_matrix(Matrix::t(peaks_near_genes),
                                 groups = colnames(peaks_near_genes),
                                 fun = "sum")
  # Keep only the genes that are in our scRNA-seq dataset
  peaks2genes <- peaks2genes[rownames(peaks2genes)
                 %in% rownames(hummus_object@assays[[gene_assay]]), ]
  # Remove rows/cols with only zeros
  peaks2genes <- peaks2genes[Matrix::rowSums(peaks2genes) != 0,
                             Matrix::colSums(peaks2genes) != 0]
  # peak-gene links
  peaks2genes <- expand.grid(rownames(peaks2genes),
                          colnames(peaks2genes))[as.vector(peaks2genes > 0), ]
  colnames(peaks2genes) <- c("gene", "peak") # set column names


  # Save peak-gene links
  store_network(network = peaks2genes,
                store_network = store_network,
                output_file = output_file,
                verbose = 1)

  # Set default names for the networks if not provided
  if (is.null(gene_multiplex_name)) {
    gene_multiplex_name <- gene_assay
  }
  if (is.null(peak_multiplex_name)) {
    peak_multiplex_name <- peak_assay
  }

  # Return atac-rna bipartite
  hummus_object@multilayer@bipartites[[bipartite_name]] <- new("bipartite",
                           "network" = peaks2genes,
                           "multiplex_left" = gene_multiplex_name,
                           "multiplex_right" = peak_multiplex_name)
  return(hummus_object)
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
  sep_Peak1 = "-",
  sep_Peak2 = "-",
  sep_RegEl1 = "-",
  sep_RegEl2 = "-"
  ) {
  # Make sure Peaks and RegEl are unique
  Peaks <- unique(Peaks)
  RegEl <- unique(RegEl)

  # convert genomic corrdinate string to GRanges object
  Peak_GRangesObj <- Signac::StringToGRanges(Peaks, 
                                             sep = c(sep_Peak1, sep_Peak2))
  RegEl_GRangesObj <- Signac::StringToGRanges(RegEl,
                                              sep = c(sep_RegEl1, sep_RegEl2))

  # find overlap between peaks and regulatory elements
  PeakOverlaps <- IRanges::findOverlaps(query = RegEl_GRangesObj,
                                        subject = Peak_GRangesObj)

  # return peaks that overlapped with regulatory element
  return(Peaks[unique(as.matrix(PeakOverlaps)[, 2])])
}