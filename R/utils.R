#' @title Motifs database class
#' MotifsDatabase object stores motifs(PFM matrices)
#' and tf2motifs (TF to motifs names mapping) data.
#'
#' @slot motifs (TFBSTools::PWMatrixList) - PFM matrices.
#' @slot tf2motifs (data.frame) - TF to motif names mapping. Columns: motif, tf.
#'
#' @name  motifs_db-class
#' @rdname motifs_db-class
#' @exportClass motifs_db
motifs_db = setClass("motifs_db",
                           representation(
                             motifs = "PWMatrixList",
                             tf2motifs = "data.frame"
                           ))

#' @title Extract TF names from scRNA data and tf2motifs
#'
#' @param species (character) - Species name. Default: "human".
#' @param l_expressed_genes (vector(character)) - List of expressed genes.
#' @param output_file (character) - Path to output file.
#' @param tf2motifs (data.frame) - TF to motifs names mapping. Columns: motif, tf.
#' @param verbose (integer) - Verbosity level. Default: 1.
#'
#' @return TFs (vector(character)) - List of TFs expressed with motifs in the database.
#' @export
#'
getTFs <- function(species="human", l_expressed_genes, output_file, tf2motifs, verbose=1){
  
  TFs <- intersect(unique(as.character(tf2motifs$tf)), l_expressed_genes)
  
  if (verbose>0){
    print(paste(length(TFs), "TFs expressed"))
  }
  
  write.table(TFs, output_file,                                                                    # Store TFs
              col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  return(TFs)
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
