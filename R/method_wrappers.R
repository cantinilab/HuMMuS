# Description: This file contains the wrapper functions for the methods that
# are used to compute the different layers of the multilayer network. The
# functions are called from the compute_*_network functions in layers.R
# For now, only the compute_atac_peak_network function has wrapper functions
# for the different methods. The other methods are still directly implemented
# in the compute_*_network functions in layers.R

#' @title Cicero wrapper function for the compute_atac_peak_network function
#'
#' @description This function is a wrapper for the compute_atac_peak_network
#' function in layers.R. It computes the peak network from scATAC-seq data
#' using Cicero. It returns a data frame with the peak network. The data frame
#' also contains the coaccess score for each edge. The coaccess score is the
#' probability that two peaks are accessible in the same cell. The coaccess 
#' score is computed by Cicero. Edges are filtered based on the coaccess score.
#' Only edges with a coaccess score > threshold are kept.
#'
#' @param hummus A hummus object
#' @param atac_assay The name of the assay containing the scATAC-seq data
#' @param genome The genome object
#' @param window The window size used by Cicero to compute the coaccess score
#' @param number_cells_per_clusters The number of cells per cluster used by
#' Cicero to compute the coaccess score
#' @param sample_num The number of samples used by Cicero to compute the
#' coaccess score
#' @param seed The seed used by Cicero to compute the coaccess score
#' @param verbose The verbosity level
#' @param threshold The threshold used to filter edges based on the coaccess
#' score
#' @param reduction_method The method used by monocle3 to reduce the dimension
#' of the scATAC-seq data before defining the pseudocells. The default is UMAP.
#'
#' @return A data frame containing the peak network
#' @export
#'
run_cicero_wrapper <- function(
    hummus,
    atac_assay,
    genome,
    window,
    number_cells_per_clusters,
    sample_num,
    seed,
    verbose,
    threshold,
    reduction_method = "UMAP"
    ) {

    # functions that need to be renamed
    int_elementMetadata <- SingleCellExperiment::int_elementMetadata
    counts <- SingleCellExperiment::counts

    # obtain chromosome sizes
    chromosome_sizes <- data.frame(V1 = genome@seqinfo@seqnames,
                                   V2 = genome@seqinfo@seqlengths)

    # Get scATAC-seq data
    scATAC <- as.matrix(hummus@assays[[atac_assay]]@counts)
    # Matrix to edgelist
    acc <- reshape2::melt(scATAC)
    colnames(acc) <- c("V1", "V2", "V3")

    # Prepare cicero input
    input_cds <- cicero::make_atac_cds(acc, binarize = TRUE) # Create CDS object
    set.seed(seed)
    # It is required that there is no empty cell
    if (length(which(colSums(as.matrix(monocle3::exprs(input_cds))) == 0)) == 0
    ) {
    # Calculating size factors using default method = mean-geometric-mean-total
      input_cds <- monocle3::estimate_size_factors(input_cds)
      # Preprocessing using LSI
      input_cds <- monocle3::preprocess_cds(input_cds, method = "LSI")
      # Dimensionality reduction using UMAP
      input_cds <- monocle3::reduce_dimension(
                                    input_cds,
                                    reduction_method = reduction_method,
                                    preprocess_method = "LSI")
    } else {
      print("Error: there is at least one cell with no signal.")
    }
    # Get reduced (UMAP) coordinates
    umap_coords <- SingleCellExperiment::reducedDims(input_cds)$UMAP
    # Compute pseudocells
    cicero_cds <- cicero::make_cicero_cds(
      input_cds,  # Create a Cicero CDS object
      reduced_coordinates = umap_coords,
      k = number_cells_per_clusters,  #number neighbors/ Default = 50
      summary_stats = NULL,         # Default
      size_factor_normalize = TRUE, # Default
      silent = FALSE)               # Default

    cicero <- cicero::run_cicero(
      cds = cicero_cds, # Infer peak-links
      genomic_coords = chromosome_sizes,
      window = window,             # Default = 5e+05
      silent = FALSE,             # Default = FALSE
      sample_num = sample_num) # Default = 100

    # Remove NAs, double edges, and edges with coaccess score <=0
    # Check for coaccess = NA
    if (length(which(is.na(cicero$coaccess))) > threshold) {
      cicero <- cicero[which(!is.na(cicero$coaccess)), ]  # Remove NAs
    }
    cicero$temp <- NA  # Helper column to check and remove double edges
    my_cols <- which(as.character(cicero$Peak1) <= as.character(cicero$Peak2))
    cicero$temp[my_cols] <- paste(cicero$Peak1[my_cols],
                                  cicero$Peak2[my_cols],
                                  sep = ";")

    my_cols <- which(as.character(cicero$Peak1) > as.character(cicero$Peak2))
    cicero$temp[my_cols] <- paste(cicero$Peak2[my_cols],
                                  cicero$Peak1[my_cols],
                                  sep = ";")

    # Sort table according to temp-column (each entry appears double)
    cicero <- cicero[with(cicero, order(temp, decreasing = TRUE)), ]
    rownames(cicero) <- c(1:dim(cicero)[1])
    A <- as.character(cicero$Peak1[seq(1, dim(cicero)[1], 2)])
    Anum <- round(cicero$coaccess[seq(1, dim(cicero)[1], 2)], 10)
    B <- as.character(cicero$Peak2[seq(2, dim(cicero)[1], 2)])
    Bnum <- round(cicero$coaccess[seq(2, dim(cicero)[1], 2)], 10)
    # length(which(A==B & Anum==Bnum)) 
    # Each edge appears twice with same coaccess score (rounded to 10 digits)
    cicero <- cicero[seq(1, dim(cicero)[1], 2), ] # Remove double edges
    cicero$temp <- NULL # Remove helper column
    cicero <- cicero[with(cicero, order(cicero$coaccess,
                                        decreasing = TRUE)), ]  # Sort
    rownames(cicero) <- c(1:dim(cicero)[1])
    cicero$Peak1 <- gsub("_", "-", cicero$Peak1)
    # Peak names 2x"-" to match bipartites
    cicero$Peak2 <- gsub("_", "-", cicero$Peak2)
    # Peak names 2x"-" to match bipartites ? 2x"-" or 2x"_"

    peak_network <- cicero[which(cicero$coaccess > threshold), ]
    # Remove edges with coaccess score <= threshold

    if (verbose > 0) {
      cat("\n", dim(peak_network)[1], "peak edges with a coaccess score >",
          threshold, "were found.\n")
    }

    # Return peak network including edges with positive coaccess score
    return(peak_network)
}


#' @title Omnipath wrapper function for the compute_tf_network function
#' @description This function is a wrapper for the compute_tf_network function
#' in layers.R. It computes the TF network from using Omnipath database.
#' It returns a data frame with the TF network. The data frame is not weighted
#' and does not contain scores for the edges.
#' @param hummus A hummus object
  # Get TF-TF interactions from Omnipath
run_omnipath_wrapper <- function(
  hummus = hummus,
  organism = organism,
  tfs = tfs,
  gene_assay = gene_assay,
  source_target = source_target,
  verbose = 1) {

  TF_PPI <- OmnipathR::import_post_translational_interactions(
    organism = organism, partners = tfs, source_target = source_target
  )

  if (verbose > 0) {
    cat("\tNumber of edges from Omnipath:", nrow(TF_PPI),
    "\nWill now be filtered to only those corresponding to specified tfs")
  }
  
  if (is.na(tfs)) {
    # Get tfs list
    tfs <- get_tfs(hummus = hummus,
              assay = gene_assay,
              store_tfs = FALSE,
              output_file = NULL,
              verbose = verbose)
  } else if (typeof(tfs) != "character") {
      stop("'tfs' argument needs to be a vector of characters
      (e.g.: c('MYC', 'JAK1')).")
  }

  # add filtering if element is not a TF expressed in the dataset
  if (source_target == "AND") {
    TF_PPI <- TF_PPI[which(TF_PPI$source_genesymbol %in% tfs &
                           TF_PPI$target_genesymbol %in% tfs), ]
  } else if (source_target == "OR") {
    TF_PPI <- TF_PPI[which(TF_PPI$source_genesymbol %in% tfs |
                           TF_PPI$target_genesymbol %in% tfs), ]
  }
  # Get only source and target columns
  tf_network <- TF_PPI[, c(3, 4)]

  if (verbose > 0) {
    cat("\tTF network construction time:", Sys.time() - a, "\n")
  }
  # Convert to data.frame from tibble
  tf_network <- as.data.frame(tf_network)

  # Check if there is any TF-TF edges otherwise add a fake node
  # and connect all TFs to it (to allow HuMMuS to run without impacting result)
  if (nrow(tf_network) == 0) {
    cat("No TF-TF edges from Omnipath for the given parameters.
        You can try to change the source_target parameter to 'OR' to get
        TF-other protein interactions. Or try to import a network  
        computed externally. Right now, a network with all TFs connected
        to a fake node is created, for HuMMuS analysis.\n It has no biological
        meaning but will allow to run the pipeline as if no edges were present.
        \n")
    run_tf_null_wrapper()

  }
  return(tf_network)
}

#' @title tf_null wrapper function for the tf_network function
#' @description This function is a wrapper for the tf_network function
#' 
#' @param hummus A hummus object
#' 
#' 
run_tf_null_wrapper <- function(
  hummus = hummus,
  organism = organism,
  tfs = tfs,
  gene_assay = gene_assay,
  verbose = 1) {
    
  if (verbose > 0) {
    cat("Creating a fake TF network with all TFs connected to a fake node.\n")
  }

  if (is.na(tfs)) {
    # Get tfs list
    tfs <- get_tfs(hummus = hummus,
              assay = gene_assay,
              store_tfs = FALSE,
              output_file = NULL,
              verbose = verbose)
  } else if (typeof(tfs) != "character") {
      stop("'tfs' argument needs to be a vector of characters
      (e.g.: c('MYC', 'JAK1')).")
  }


  FAKE_NODE <- "fake_node"
  tf_network <- data.frame(colnames = c("source", "target"))
  for (tf in tfs) {
    tf_network <- rbind(tf_network, data.frame(source = tf, target = FAKE_NODE))
  }
  return(tf_network)
}