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
get_tfs <- function(
    hummus,
    assay = NULL,
    store_tfs = TRUE,
    output_file = NULL,
    verbose = 0
    ) {
  # Check if the hummus object has motifs_db slot
  if (is.null(hummus@motifs_db)) {
    stop("The hummus object does not have a motifs_db slot")
  }
  
  # Check if the assay is present in the seurat object
  if (! is.null(assay)) {
    if (!assay %in% names(hummus@assays)) {
        stop("The gene assay is not present in the seurat object")
    }
    # Get the expressed genes
    expr_genes <- rownames(hummus@assays[[assay]])
    tfs <- intersect(unique(as.character(hummus@motifs_db@tf2motifs$tf)),
                                        expr_genes)
    if (verbose > 0) {
        cat("\t", length(tfs), "TFs expressed\n")
        }
  } else { # If no assay is provided, get all TFs with motifs
    tfs <- unique(as.character(hummus@motifs_db@tf2motifs$tf))
    if (verbose > 0) {
      cat("\t", length(tfs), "TFs with motif. No check if expressed or not.\n")
    }
  }
  # Store TFs in a file if specified
  if (store_tfs) {
    if (is.null(output_file)) {
      stop("Please provide an output file name")
    }
    write.table(tfs, output_file, # Store TFs
              col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  }

  return(tfs)
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
) {
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
        x <- gsub(pattern = ', ', replacement='_', x=x, fixed=TRUE)
        x <- gsub(pattern = '_sep = \"_\")',
                  replacement = "",
                  x = x,
                  fixed = TRUE)
        return(x)
    })

    result <- result * values
    if(isTRUE(response > 0))
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

    if ('character'%in%class(fun)) {
        fun <- summary_fun[[fun]]
    }

    if (length(groups) == nrow(x)) {
        agg_mat <- sapply(levels(factor(groups)), function(g){
            chunk <- x[which(groups==g), ]
            if (is.null(dim(chunk))){
                return(chunk)
            } else {
                return(fun(chunk))
            }
        })
        agg_mat <- Matrix::Matrix(agg_mat, sparse=TRUE)
    } else if (length(groups) <= 1) {
        agg_mat <- fun(x)
        agg_mat <- Matrix::Matrix(agg_mat, sparse=TRUE)
        colnames(agg_mat) <- groups
        rownames(agg_mat) <- colnames(x)
    } else {
        stop('Length of groups must be either nrow(x) or 1.')
    }
    return(Matrix::t(agg_mat))
}