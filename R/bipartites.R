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
Bipartite_TFs2Peaks <- function(
  tfs,
  peaks,
  peak_sep1 = ":",
  peak_sep2 = "-",
  genome,
  gene.range,
  motifs,
  tf2motifs,
  store_bipartite = TRUE,
  output_file = None,
  verbose = 1) {

  # Build up object to determine TF-peak links and peak-gene links
  rna <- data.frame(features=tfs) # List of genes present in scRNA
  rna[, c("dummy_cell1", "dummy_cell2")] <- 1
  rownames(rna) <- rna$features
  rna$features <- NULL
  atac  <- data.frame(features=peaks) # List of peaks present in scATAC
  atac[, c("dummy_cell1", "dummy_cell2")] <- 1
  rownames(atac) <- atac$features
  atac$features <- NULL
  seurat <- CreateSeuratObject(rna)
  # Create seurat object with genes
  seurat[["peaks"]] <- CreateChromatinAssay(atac, sep = c(peak_sep1,peak_sep2))
  # Combine genes and peaks in a seurat object
  Annotation(seurat@assays$peaks) <- gene.range
  # Add genome annotations to seurat object

  motif_pos <- Signac::AddMotifs(
    object = seurat[["peaks"]],
    genome = genome,
    pfm = motifs,  #add verbose options
  )

  ## The 17 following lines are inspired from the Pando package :
  # https://github.com/quadbiolab/Pando/blob/main/R/regions.R
  # Add TF info for motifs
  if (verbose > 0) {
    print("Adding TF info", verbose = verbose)
  }
  if (!is.null(tf2motifs)) {
    tf2motifs <- dplyr::tibble(tf2motifs) #if error check other than dplyr::
  } else {
    utils::data(tf2motifs, envir = environment())
  }

  # Spread dataframe to sparse matrix
  tf2motifs <- tf2motifs %>%
    dplyr::select("motif" = 1, "tf" = 2) %>%
    dplyr::distinct() %>%
    dplyr::mutate(val = 1) %>%
    tidyr::pivot_wider(names_from = "tf", values_from = val, values_fill = 0) %>%
    tibble::column_to_rownames("motif") %>%
    as.matrix() %>%
    Matrix::Matrix(sparse = TRUE)
  tfs_use <- intersect(rownames(GetAssay(seurat, "RNA")), colnames(tf2motifs))
  if (length(tfs_use) == 0) { # If no TFs are found in the dataset
    stop("None of the provided TFs were found in the dataset.
    Consider providing a custom motif-to-TF map as `motif_tfs`")
  }

  TFs_Peaks <- motif_pos@motifs@data %*% tf2motifs[, tfs_use] # Get TF peak links
  TFs_Peaks <- TFs_Peaks[, colnames(TFs_Peaks) %in% tfs]
   # Keep only the TFs that are in our scRNA-seq dataset
  tfs2peaks <- expand.grid(rownames(TFs_Peaks),
                           colnames(TFs_Peaks))[as.vector(TFs_Peaks > 0), ]
                          # TF-peak links
  colnames(tfs2peaks) <- c("peak", "TF")     # set column names

  if (store_bipartite) { # Save TF-peak links
    write.table(tfs2peaks,
                output_file,
                col.names = TRUE,
                row.names = FALSE,
                quote = FALSE,
                sep = "\t")
  }

  return(tfs2peaks) # Return TF-peak links
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
Bipartite_Peaks2Genes <- function(
  seurat_object,
  gene_assay = 'RNA',
  peak_assay = 'peaks',
  peak_to_gene_method = "Signac",
  upstream = 500,
  downstream = 500,
  only_tss = TRUE,
  store_bipartite = TRUE,
  output_file = None) {
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
    write.table(peaks2genes,
                output_file,
                col.names = TRUE,
                row.names = FALSE,
                quote = FALSE,
                sep = "\t")
  }

  # Return list the two edgelists containing TF-peak and peak-gene links
  return(peaks2genes)
}
