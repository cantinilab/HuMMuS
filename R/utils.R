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
motifs_db <- setClass("motifs_db",
                           representation(
                             motifs = "PWMatrixList",
                             tf2motifs = "data.frame",
                             tfs = "NULL"
                           ))
setMethod("show", "motifs_db",
  function(object) {
    cat(
      paste("Motifs database object with :\n- ",
          length(object@motifs), "motifs\n- ",
          length(unique(object@tf2motifs$tf)), " TFs\n- ",
          nrow(object@tf2motifs), "TF to motif names mapping"
          )
      )
  })

multiplex <- setClass(Class = "multiplex",
                       slots = c(
                         "networks" = "list", # List of networks
                         "features" = "vector"# Vector of features
                        )
                      )
setMethod("show", "multiplex",
  function(object) {
    cat(
      paste("Multiplex of ", length(object@networks),
      " networks with", length(object@features), "features.\n",
      "Networks names: ", paste(names(object@networks), collapse = ", "))
      )
  })

bipartite <- setClass(Class = "bipartite",
                       slots = c(
                      "network" = "data.frame", # Bipartite network (edge list)
                      "multiplex_left" = "character", # left features' multiplex
                      "multiplex_right" = "character" # right features multiplex
                        )
                      )
setMethod("show", "bipartite",
  function(object) {
    cat(
      paste("Bipartite network of ", nrow(object@network), " edges.\n",
      "Multiplexes names: ", object@multiplex_left,
      " and ", object@multiplex_right, "\n")
      )
  })

multilayer <- setClass(Class = "multilayer",
                       slots = c(
                        "bipartites" = "list", # Bipartite networks
                        "multiplex" = "list", # Multiplex networks
                        "config" = "list" # Parameters to create the hmln
                        )                 # representation of a yaml file
                      )
setMethod("show", "multilayer",
  function(object) {
    cat(
      paste("Multilayer network of ",
      length(object@bipartites), " bipartite networks and ",
      length(object@multiplex), " multiplex networks.\n",
      "\n- Multiplex networks names: ", paste(names(object@multiplex),
                                          collapse = ", "),
      "\n- Bipartite networks names: ", paste(names(object@bipartites),
                                          collapse = ", "), "\n"
      )
    )
  })

#' The hummus_object class
#'
#' The SeuratPlus object is an extended \code{Seurat} object
#' for the storage and analysis of a heterogeneous mutlilayer network
#'
#' @slot multilayer
#'
#' @name hummus_object-class
#' @rdname hummus_object-class
#' @exportClass hummus_object
#' @concept assay
hummus_object <- setClass(
    Class = "hummus_object",
    contains = "Seurat",
    slots = list(
        "multilayer" = "multilayer",
        "motifs_db" = "motifs_db"
    )
)
setMethod("show", "hummus_object",
  function(object) {
    object <- SeuratObject::UpdateSlots(object = object)
    assays <- SeuratObject::FilterObjects(object = object,
                                          classes.keep = "Assay")
    nfeatures <- sum(vapply(
      X = assays,
      FUN = function(x) {
        return(nrow(x = object[[x]]))
      },
      FUN.VALUE = integer(length = 1L)
    ))
    num.assays <- length(x = assays)

    cat(
      paste("Hummus object containing a multilayer object :",
      "\n- Multiplex networks names: ",
      paste(names(object@multilayer@multiplex),
                  ollapse = ", "),
      "\n- Bipartite networks names: ",
      paste(names(object@multilayer@bipartites),
                  collapse = ", "), "\n")
    )
    cat(
      nfeatures,
      "features across",
      ncol(x = object),
      "samples within",
      num.assays,
      ifelse(test = num.assays == 1, yes = "assay", no = "assays"),
      "\n"
    )
    cat(
      "Active assay:",
      SeuratObject::DefaultAssay(object = object),
      paste0('(', nrow(x = object), " features, ",
      length(x = SeuratObject::VariableFeatures(object = object)), " variable features)")
    )
    other.assays <- assays[assays != SeuratObject::DefaultAssay(object = object)]
    if (length(x = other.assays) > 0) {
      cat(
        '\n',
        length(x = other.assays),
        'other',
        ifelse(test = length(x = other.assays) == 1, yes = 'assay', no = 'assays'),
        'present:',
        strwrap(x = paste(other.assays, collapse = ', '))
      )
    }
    reductions <- SeuratObject::FilterObjects(object = object, classes.keep = 'DimReduc')
    if (length(x = reductions) > 0) {
      cat(
        '\n',
        length(x = reductions),
        'dimensional',
        ifelse(test = length(x = reductions) == 1, yes = 'reduction', no = 'reductions'),
        'calculated:',
        strwrap(x = paste(reductions, collapse = ', '))
      )
    }
    fovs <- SeuratObject::FilterObjects(object = object, classes.keep = 'FOV')
    if (length(x = fovs)) {
      cat(
        '\n',
        length(x = fovs),
        'spatial',
        ifelse(test = length(x = fovs) == 1L, yes = 'field', no = 'fields'),
        'of view present:',
        strwrap(x = paste(fovs, sep = ', '))
      )
    }
    images <- SeuratObject::FilterObjects(object = object, classes.keep = 'SpatialImage')
    images <- setdiff(x = images, y = fovs)
    if (length(x = images)) {
      cat(
        '\n',
        length(x = images),
        ifelse(test = length(x = images) == 1L, yes = 'image', no = 'images'),
        'present:',
        strwrap(x = paste(images, collapse = ', '))
      )
    }
    cat('\n')
  }
)


#' @title Extract TF names from scRNA data and tf2motifs
#'
#' @param species (character) - Species name. Default: "human".
#' @param genes (vector(character)) - List of expressed genes.
#' @param output_file (character) - Path to output file.
#' @param tf2motifs (data.frame) - TF to motifs names mapping.
#' Columns: motif, tf.
#' @param verbose (integer) - Verbosity level. Default: 1.
#'
#' @return TFs (vector(character)) - List of TFs expressed with motifs.
#' @export
#'
get_tfs <- function(hummus,
                    assay,
                    store_tfs = TRUE,
                    output_file = NULL,
                    verbose = 1) {
  # Check if the assay is present in the seurat object
  if (!assay %in% names(hummus@assays)) {
    stop("The gene assay is not present in the seurat object")
  }

  # Check if the hummsu object has motifs_db 
  else if (is.null(hummus@motifs_db)) {
    stop("The hummus object does not have a motifs_db slot")
  }

  expr_genes = rownames(hummus@assays[[assay]])
  tfs <- intersect(unique(as.character(hummus@motifs_db@tf2motifs$tf)),
                                       expr_genes)
  if (verbose > 0) {
    print(paste(length(tfs), "TFs expressed"))
  }

  if (store_tfs) {
    if (is.null(output_file)) {
      stop("Please provide an output file name")
    }
    write.table(tfs, output_file, # Store TFs
              col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  }

  return(tfs)
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



# Code from Pando github.com/quadbiolab/Pando
#' @import sparseMatrixStats
summary_fun <- list(
    'mean' = sparseMatrixStats::colMeans2,
    'median' = sparseMatrixStats::colMedians,
    'max' = sparseMatrixStats::colMaxs,
    'min' = sparseMatrixStats::colMins,
    'count' = sparseMatrixStats::colCounts,
    'any' = sparseMatrixStats::colAnys,
    'all' = sparseMatrixStats::colAlls,
    'sd' = sparseMatrixStats::colSds,
    'mad' = sparseMatrixStats::colMads
)

#' Copy of the aggregate.Matrix function from the Matrix.utils package,
#' since this is off CRAN and does not seem to be maintained anymore
#' @keyword internal
#'
fast_aggregate <- function(
    x,
    groupings = NULL,
    form = NULL,
    fun = 'sum',
    ...
){
    if (!is(x,'Matrix')){
        x <- Matrix(as.matrix(x), sparse=TRUE)
    }
    if (fun=='count'){
        x <- x!=0
    }
    groupings2 <- groupings
    if (!is(groupings2, 'data.frame')){
        groupings2 <- as.data.frame(groupings2)
    }
    groupings2 <- data.frame(lapply(groupings2, as.factor))
    groupings2 <- data.frame(interaction(groupings2, sep='_'))
    colnames(groupings2) <- 'A'
    if (is.null(form)){
        form <- as.formula('~0+.')
    }
    form <- as.formula(form)
    mapping <- dMcast(groupings2, form)
    colnames(mapping) <- substring(colnames(mapping), 2)
    result <- Matrix::t(mapping) %*% x
    if (fun=='mean'){
        result@x <- result@x/(fast_aggregate(x, groupings2, fun='count'))@x
    }
    attr(result,'crosswalk') <- grr::extract(groupings, match(rownames(result), groupings2$A))
    return(result)
}

#' Copy of the dMcast function from the Matrix.utils package,
#' since this is off CRAN and does not seem to be maintained anymore
#' @keyword internal
#'
dMcast <- function(
    data,
    formula,
    fun.aggregate = 'sum',
    value.var = NULL,
    as.factors = FALSE,
    factor.nas = TRUE,
    drop.unused.levels = TRUE
){
    values <- 1
    if (!is.null(value.var)){
        values <- data[,value.var]
    }
    alltms <- terms(formula, data=data)
    response <- rownames(attr(alltms, 'factors'))[attr(alltms, 'response')]
    tm <- attr(alltms, "term.labels")
    interactionsIndex <- grep(':', tm)
    interactions <- tm[interactionsIndex]
    simple <- setdiff(tm, interactions)
    i2 <- strsplit(interactions,':')
    newterms <- unlist(lapply(i2, function (x) paste("paste(", paste(x, collapse=','), ",", "sep='_'",")")))
    newterms <- c(simple, newterms)
    newformula <- as.formula(paste('~0+', paste(newterms, collapse='+')))
    allvars <- all.vars(alltms)
    data <- data[, c(allvars), drop=FALSE]
    if (as.factors)
        data <- data.frame(lapply(data, as.factor))
    characters <- unlist(lapply(data, is.character))
    data[,characters] <- lapply(data[, characters,drop=FALSE], as.factor)
    factors <- unlist(lapply(data, is.factor))
    # Prevents errors with 1 or fewer distinct levels
    data[,factors] <- lapply(data[,factors,drop=FALSE],function (x)
    {
        if (factor.nas){
            if (any(is.na(x))){
                levels(x) <- c(levels(x),'NA')
                x[is.na(x)] <- 'NA'
            }
        }
        if (drop.unused.levels){
            if (nlevels(x)!=length(na.omit(unique(x)))){
                x <- factor(as.character(x))
            }
        }
        y <- contrasts(x, contrasts=FALSE, sparse=TRUE)
        attr(x, 'contrasts') <- y
        return(x)
    })
    # Allows NAs to pass
    attr(data,'na.action') <- na.pass
    result <- Matrix::sparse.model.matrix(newformula, data,drop.unused.levels = FALSE, row.names=FALSE)
    brokenNames <- grep('paste(', colnames(result), fixed = TRUE)
    colnames(result)[brokenNames] <- lapply(colnames(result)[brokenNames], function (x) {
        x <- gsub('paste(', replacement='', x=x, fixed = TRUE)
        x <- gsub(pattern=', ', replacement='_', x=x, fixed=TRUE)
        x <- gsub(pattern='_sep = \"_\")', replacement='', x=x, fixed=TRUE)
        return(x)
    })

    result <- result*values
    if(isTRUE(response>0))
    {
        responses=all.vars(terms(as.formula(paste(response,'~0'))))
        result <- fast_aggregate(result, data[, responses,drop=FALSE], fun=fun.aggregate)
    }
    return(result)
}


#' Aggregate matrix over groups
#'
#' @import sparseMatrixStats
#'
#' @param groups A character vector with the groups to aggregate over.
#' @param fun The summary function to be applied to each group.
#'
#' @return A summary matrix.
#'
#' @export
aggregate_matrix <- function(
    x,
    groups = NULL,
    fun = 'mean'
){
    if (length(groups) == nrow(x) & 'character'%in%class(fun)){
        if (fun%in%c('count', 'sum')){
            agg_mat <- fast_aggregate(x=x, groupings=groups, fun=fun)
            return(agg_mat)
        }

        if (fun=='mean'){
            group_counts <- as.numeric(table(groups))
            agg_mat <- fast_aggregate(x=x, groupings=groups, fun='sum')
            agg_mat <- agg_mat / group_counts
            return(agg_mat)
        }
    }

    if ('character'%in%class(fun)){
        fun <- summary_fun[[fun]]
    }

    if (length(groups) == nrow(x)){
        agg_mat <- sapply(levels(factor(groups)), function(g){
            chunk <- x[which(groups==g), ]
            if (is.null(dim(chunk))){
                return(chunk)
            } else {
                return(fun(chunk))
            }
        })
        agg_mat <- Matrix::Matrix(agg_mat, sparse=TRUE)
    } else if (length(groups) <= 1){
        agg_mat <- fun(x)
        agg_mat <- Matrix::Matrix(agg_mat, sparse=TRUE)
        colnames(agg_mat) <- groups
        rownames(agg_mat) <- colnames(x)
    } else {
        stop('Length of groups must be either nrow(x) or 1.')
    }
    return(Matrix::t(agg_mat))
}
