#' Compute TF network
#'
#' Compute a protein-protein interaction layer from Omnipath request that will represent tf cooperativity.
#' This network is the top-layer of HuMMuS multilayer.
#' 
#' @param organism (integer)  - Specie identifier from Omnipath to fetch specific interactions
#' @param tfs vector(character) - List of tfs considered.
#' @param store_network (bool) - Save the network directly (\code{TRUE},
#'  default) or return without saving on disk (\code{FALSE}).
#' @param output_file (character) - Name of the output_file (if store_network == \code{TRUE}).
#' @param source_target ('AND'|'OR') - Fetch only the interactions involving
#'  two considered tfs (\code{'AND'}), or one considered tfs and any other element (\code{'OR'})
#'
#' @return (data.frame) - Return list of network interactions between tfs (or larger set of protein if source_target==\code{'OR'})
#' @export
#'
#' @examples TO DO. Same than UNIT test.
compute_tf_network <- function(hummus = null, # Hummus object
                               organism = 9606, # Human by default
                               tfs = NA, # List of tfs considered.
                               gene_assay = NULL, # Name of the assay to get tfs from
                                                  # if tfs is not provided
                               store_network = FALSE, # Save the network on disk (TRUE, default)
                               output_file = NULL, # Name of the output_file (if store_network == TRUE)
                               source_target = "AND", # 'AND' | 'OR'
                               only_expressed = TRUE, # TRUE | FALSE
                               multiplex_name = NULL, # Name of the multiplex to add the network to
                               tf_network_name = "TF_network", # Name of the network in the multiplex
                               verbose = 1) {

  if (verbose > 0) {
    cat("Computing TF network...\n")
    a <- Sys.time()
  }
  TF_PPI <- OmnipathR::import_post_translational_interactions(
    organism = organism, partners = tfs, source_target = source_target
  )

  if (!is.null(gene_assay)) {
    tfs <- get_tfs(hummus = hummus,
            assay = gene_assay,
            store_tfs = FALSE,
            output_file = NULL,
            verbose = verbose)
  }
  # add filtering if element is not a TF expressed in the dataset
  if (source_target == "AND") {
    TF_PPI <- TF_PPI[which(TF_PPI$source %in% tfs &
                           TF_PPI$target %in% tfs), ]
  } else if (source_target == "OR") {
    TF_PPI <- TF_PPI[which(TF_PPI$source %in% tfs |
                           TF_PPI$target %in% tfs), ]
  }

  Network_TF <- TF_PPI[, c(3, 4)]

  if (verbose > 0) {
    cat("\tTF network construction time:", Sys.time() - a, "\n")
  }

  if (store_network == TRUE) {
    if (verbose > 0) {
      cat("\tStoring TF network in file : ", output_file, "\n")
    }
    if (is.null(output_file)) {
      stop("\tPlease provide an output file name if you want to store the network.\n")
    }
    write.table(Network_TF, output_file, # Store edgelist
                col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  }

  if (is.null(hummus)) {
    cat('\tno hummus object provided\n')
    return(Network_TF)
  } else if (is.null(hummus@multilayer@multiplex[[multiplex_name]])) {
    if (verbose > 0) {
      cat("\tCreating new multiplex : ", multiplex_name, "\n")
    }
    hummus@multilayer@multiplex[[multiplex_name]] <- new("multiplex")
    hummus@multilayer@multiplex[[multiplex_name]] <- add_network(
                    hummus@multilayer@multiplex[[multiplex_name]],
                    network = Network_TF,
                    network_name = tf_network_name)
  } else {
    if (verbose > 0) {
      cat("\tAdding network to existing multiplex : ", multiplex_name, "\n")
    }
    hummus <- add_network(hummus,
                          multiplex_name = multiplex_name,
                          network = Network_TF,
                          network_name = tf_network_name)
  }

  return(hummus)
}


#' Compute gene netwok from scRNA-seq data
#'
#' This function will create a network from rna data (or in theory any data
#' wtih genes as features).
#' Different method should be implemented at some point (any suggestion is welcomed ! :) ),
#' for now Genie3 is still the reference and only method available
#'
#' Method descriptions :
#'  1. Genie3
#'      Use tree random forest to infer regulatory networks :
#'      https://bioconductor.org/packages/release/bioc/html/GENIE3.html
#'
#' @param scRNA (data.frame) - Expression matrix (cells*genes)
#' @param tfs vector(character) - List of tfs considered.
#' @param method (character) - Method used to infer network edges.
#' * \code{'Genie3'} - TO DO.
#' @param store_network (bool) - Save the network directly (\code{TRUE},
#'  default) or return without saving on disk (\code{FALSE}).
#' @param output_file (character) - Name of the output_file (if store_network == \code{TRUE}).
#' @param threshold (interger, default 0) - Minimal threshold to select tf-gene edges.
#' @param number_cores (interger, default 1) - Number of thread that should be used for the parallelizable methods.
#' @param verbose (integer) - Display function messages. Set to 0 for no message displayed, >= 1 for more details.
#'
#' @return (data.frame) - Return list of network interactions between genes
#' @export
#'
#' @examples TO DO. Same than UNIT test.
compute_gene_network <- function(hummus,
                                 gene_assay = "RNA",
                                 tfs = NULL,
                                 method = "GENIE3",
                                 store_network = FALSE,
                                 output_file,
                                 threshold = 0.0,
                                 number_cores = 1,
                                 verbose = 1,
                                 multiplex_name = NULL) {
  if (method == "GENIE3") {
    if (verbose > 0) {
      cat("Computing gene network with ", method, " ...\n")
      a <- Sys.time()
    }
    # Get tfs list
    if (verbose > 0 && is.null(tfs)) {
      cat("\tNo TFs list provided, fetching from hummus object...\n")
    }
    tfs <- get_tfs(hummus = hummus,
            assay = gene_assay,
            store_tfs = FALSE,
            output_file = NULL,
            verbose = verbose)

    # infer network
    weightMat <- GENIE3::GENIE3(as.matrix(hummus[[gene_assay]]@counts),
                               regulators = tfs,
                               nCores = number_cores)

    if (verbose > 0) {
      cat("\tGene network construction time:", Sys.time() - a, "\n")
    }
    # get edge list
    linkList <- GENIE3::getLinkList(weightMat)
    gene_network <- linkList[which(linkList$weight > threshold), ]


    features <- unique(c(unique(gene_network$regulatoryGene),
                      unique(gene_network$targetGene)))

  # Save gene network
  if (store_network) {
    if (is.null(output_file)) {
      stop("Please provide an output file name if you want to store the network.")
    }
    if (verbose > 0) {
      cat("\tStoring gene network in file : ", output_file, "\n")
    }
    write.table(gene_network,
                output_file,
                col.names = TRUE,
                row.names = FALSE,
                quote = FALSE,
                sep = "\t")
  }

  # TODO : add other methods
  } else {
    stop(cat("Method not implemented yet, choose between GENIE3 and..",
    "that's it for now.\n but you can always compute the network",
    "independently and add it to the hummus object."))
  }

  # Check if multiplex name provided
  if (is.null(multiplex_name)) {
    multiplex_name <- gene_assay
  }
  # Create a new multiplex of the hummus multilayer if does not exist already
  if (is.null(hummus@multilayer@multiplex[[multiplex_name]])) {
    if (verbose > 0) {
      cat("\tCreating new multiplex : ", multiplex_name, "\n")
    }
    hummus@multilayer@multiplex[[multiplex_name]] <- new("multiplex")
    hummus@multilayer@multiplex[[multiplex_name]] <-
              add_network(hummus@multilayer@multiplex[[multiplex_name]],
                          network = gene_network,
                          network_name = paste0(multiplex_name, "_", method))
  } else {paste
    if (verbose > 0) {
      cat("\tAdding network to existing multiplex : ", multiplex_name, "\n")
    }
    # Add network to existing multiplex if exists already
    hummus <- add_network(hummus,
                        multiplex_name = multiplex_name,
                        network = gene_network,
                        network_name = paste0(multiplex_name, "_", method))
  }
  # Return hummus object
  return(hummus)
}

#' Compute peak network from scATAC-seq data
#'
#' This function will create a network from atac data (or in theory any data
#' wtih peaks coordinates as features).
#' Different method should be implemented at some point (e.g. RENIN),
#' for now Cicero is still the reference and only method available
#'
#' Method descriptions :
#'  1. Cicero
#'      Use patial corelation between peaks that are in a given window (e.g. :
#'      less distant than 500K base pairs)
#'
#' @param scATAC TODO
#' @param genome TODO
#' @param method TODO
#' @param store_network (bool) - Save the network directly (\code{TRUE},
#'  default) or return without saving on disk (\code{FALSE}).
#' @param output_file (character) - Name of the output_file (if store_network == \code{TRUE}).
#' @param threshold (interger, default 0) - Minimal threshold to select tf-gene edges.
#' @param number_cells_per_clusters (integer) - Number of cells grouped by territory to define pseudocells
#' @param sample_num (integer | Cicero) - TO DO.
#' @param seed ( __ ) - Random seed used by Cicero.
#' @param verbose (integer) - Display function messages. Set to 0 for no message displayed, >= 1 for more details.
#' @param window (interger) - Size of window to consider potential cis-regulatory cooperations
#'
#' @return (data.frame) - Return list of network interactions between peaks
#' @export
#'
#' @examples TO DO. Same than UNIT test.
compute_atac_peak_network <- function(
    scATAC,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    method = "cicero",
    store_network = TRUE,
    output_file,
    threshold = 0.0,
    number_cells_per_clusters = 50,
    sample_num = 100,
    seed = 2021,
    verbose = 1,
    window = 5e+05
    ) {

  if (method=="cicero") {
      if (!requireNamespace("cicero", quietly = TRUE)) {
        stop("Please install cicero.\n",
         "https://cole-trapnell-lab.github.io/cicero-release/docs_m3/")
      }
    int_elementMetadata <- SingleCellExperiment::int_elementMetadata

  # obtain chromosome sizes
    chromosome_sizes <- data.frame(V1 = genome@seqinfo@seqnames,
                                   V2 = genome@seqinfo@seqlengths)
    # Matrix to edgelist
    ACC <- reshape2::melt(as.matrix(scATAC))
    colnames(ACC) <- c("V1", "V2", "V3")
    ACC$V1 <- gsub("_", "-", ACC$V1)
    # Run cicero
    input_cds <- cicero::make_atac_cds(ACC, binarize = TRUE) # Create CDS object
    set.seed(seed)
    print(input_cds)
    # It is required that there is no empty cell
    if (length(which(colSums(as.matrix(Biobase::exprs(input_cds))) == 0)) == 0)
    {
      input_cds <- estimate_size_factors(input_cds)   # Calculating size factors using default method = mean-geometric-mean-total
      input_cds <- preprocess_cds(input_cds, method = "LSI")  # Preprocessing using LSI
      input_cds <- cicero::reduce_dimension(input_cds, reduction_method = 'UMAP',  # Dimensionality reduction
                                    preprocess_method = "LSI")
    }
    else{
      print("Error: there is at least one cell with no signal.")
    }
    print(input_cds)
    umap_coords <- reducedDims(input_cds)$UMAP # Get reduced (UMAP) coordinates
    cicero_cds <- cicero::make_cicero_cds(input_cds,  # Create a Cicero CDS object
                                  reduced_coordinates = umap_coords,
                                  k = number_cells_per_clusters,  #number neighbors         # Default = 50
                                  summary_stats = NULL,         # Default
                                  size_factor_normalize = TRUE, # Default
                                  silent = FALSE)               # Default
    print("test 2")
    a = Sys.time()
    cicero <- cicero::run_cicero(cds = cicero_cds,                                   # Infer peak-links
                         genomic_coords = chromosome_sizes, # mouse.mm10.genome,
                         window = window,             # Default
                         silent = FALSE,             # Default
                         sample_num = sample_num) # Default = 100
    print("test 3")
    if (verbose > 0) {
      cat("Peak network construction time:", Sys.time() - a)
      }
    if (store_network ){
      write.table(cicero,
                  cicero_nw_filename,
                  col.names = TRUE,
                  row.names = FALSE,
                  quote = FALSE,
                  sep = "\t")
    }
    # Remove NAs, double edges, and edges with weight <=0
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
      cat(dim(peak_network)[1], "peak edges with weight > ", as.string(threshold))}

    if (store_network) {
      write.table(peak_network,
                  output_filename,
                  col.names = FALSE,
                  row.names = FALSE,
                  quote = FALSE,
                  sep = "\t") # Store edgelist
    }

    # Return peak network including edges with positive coaccess score
    print("test")
    return(peak_network)
  } else {
    return("no peak network computed!")
  }
}