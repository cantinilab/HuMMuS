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
    ACC <- reshape2::melt(scATAC)
    colnames(ACC) <- c("V1", "V2", "V3")

    # Prepare cicero input
    input_cds <- cicero::make_atac_cds(ACC, binarize = TRUE) # Create CDS object
    set.seed(seed)
    # It is required that there is no empty cell
    if (length(which(colSums(as.matrix(monocle3::exprs(input_cds))) == 0)) == 0)
    {
      # Calculating size factors using default method = mean-geometric-mean-total
      input_cds <- monocle3::estimate_size_factors(input_cds)
      # Preprocessing using LSI
      input_cds <- monocle3::preprocess_cds(input_cds, method = "LSI")
      # Dimensionality reduction using UMAP
      input_cds <- monocle3::reduce_dimension(
                                    input_cds,
                                    reduction_method = reduction_method,
                                    preprocess_method = "LSI")
    }
    else{
      print("Error: there is at least one cell with no signal.")
    }
    # Get reduced (UMAP) coordinates
    umap_coords <- SingleCellExperiment::reducedDims(input_cds)$UMAP
    # Compute pseudocells
    cicero_cds <- cicero::make_cicero_cds(input_cds,  # Create a Cicero CDS object
                                  reduced_coordinates = umap_coords,
                                  k = number_cells_per_clusters,  #number neighbors         # Default = 50
                                  summary_stats = NULL,         # Default
                                  size_factor_normalize = TRUE, # Default
                                  silent = FALSE)               # Default
    a = Sys.time()
    cicero <- cicero::run_cicero(cds = cicero_cds,                                   # Infer peak-links
                         genomic_coords = chromosome_sizes,
                         window = window,             # Default = 5e+05
                         silent = FALSE,             # Default = FALSE
                         sample_num = sample_num) # Default = 100

    if (verbose > 0) {
      cat("Peak network construction time:", Sys.time() - a)
      }

    # Remove NAs, double edges, and edges with coaccess score <=0
    # Check for coaccess = NA
    if (length(which(is.na(cicero$coaccess))) > threshold) {
      cicero <- cicero[which(!is.na(cicero$coaccess)),]  # Remove NAs
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
    cicero <- cicero[with(cicero, order(coaccess, decreasing = TRUE)), ]  # Sort
    rownames(cicero) <- c(1:dim(cicero)[1])
    cicero$Peak1 <- gsub("_", "-", cicero$Peak1) # Peak names 2x"-" to match bipartites
    cicero$Peak2 <- gsub("_", "-", cicero$Peak2) # Peak names 2x"-" to match bipartites########? 2x"-" or 2x"_"

    peak_network <- cicero[which(cicero$coaccess > threshold), ]
    # Remove edges with coaccess score <= threshold

    if (verbose > 0) {
      cat(dim(peak_network)[1], "peak edges with a coaccess score >",
          threshold, "were found.\n")
    }

    # Return peak network including edges with positive coaccess score
    return(peak_network)
}


store_network <- function(
    network,
    store_network,
    output_file,
    verbose = 1) {

  if (store_network) {
    if (is.null(output_file)) {
      stop("Please provide an output file name",
           " if you want to store the network.")
    }
    if (verbose > 0) {
      cat("\tStoring network in file : ", output_file, "\n")
    }
    write.table(network,
                output_file,
                col.names = TRUE,
                row.names = FALSE,
                quote = FALSE,
                sep = "\t")
  }
}

add_network <- function(
  object,
  network,
  network_name,
  multiplex_name = NULL,
  directed = FALSE,
  weighted = FALSE,
  verbose = 1) {

  # Check if object is a multiplex, a multilayer or an hummus object
  if (inherits(object, "multiplex")) {
    multiplex <- object
  } else if (inherits(object, "multilayer") ) {
    # Check if multiplex_name is NULL
    if (is.null(multiplex_name)) {
      stop("You need to specify the multiplex name.")
    }
    # Check if multiplex_name already exists
    if (!(multiplex_name %in% names(object@multilayer@multiplex))) {
      if (verbose > 0) {
        cat("\tCreating new multiplex : ", multiplex_name, "\n")
      }
      # Create new multiplex if not
      object@multilayer@multiplex[[multiplex_name]] <- new("multiplex")
    }
    # Get working multiplex
    multiplex <- object@multilayer@multiplex[[multiplex_name]]
  } else if (inherits(object, "hummus_object")) {
    # Check if multiplex_name is NULL
    if (is.null(multiplex_name)) {
      stop("You need to specify the multiplex name.")
    }
    # Check if multiplex_name already exists
    if (!(multiplex_name %in% names(object@multilayer@multiplex))) {
      if (verbose > 0) {
        cat("\tCreating new multiplex : ", multiplex_name, "\n")
      }
      # Create new multiplex if not
      object@multilayer@multiplex[[multiplex_name]] <- new("multiplex")
    }
    # Get working multiplex
    multiplex <- object@multilayer@multiplex[[multiplex_name]]

  } else {
    stop("Object is not a multiplex, a multilayer nor an hummus object.")
  }

  # Check if network name already exists in the multiplex
  if (network_name %in% names(multiplex@networks)) {
    stop("Network name already exists in the multiplex.")
  }

  # Check if there is features in common
  features <- unique(c(unique(network[, 1]), unique(network[, 2])))
  if (length(intersect(features, multiplex@features)) == 0
      && length(multiplex@features) != 0) {
    stop(cat("There is no features in common.",
      "Check if there is a mistake in the features names",
      " or if you want to create a new multiplex instead."))
  }

  # Add network
  multiplex@networks[[network_name]] <- network
  multiplex@features <- unique(c(multiplex@features, features))
  multiplex@directed[[network_name]] <- directed
  multiplex@weighted[[network_name]] <- weighted

  # Return object
  if (inherits(object, "multiplex")) {
    return(multiplex)
  } else if (inherits(object, "multilayer")) {
    object@multiplex[[multiplex_name]] <- multiplex
    return(object)
  } else if (inherits(object, "hummus_object")) {
    object@multilayer@multiplex[[multiplex_name]] <- multiplex
    return(object)
  }
}