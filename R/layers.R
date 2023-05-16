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
  tf_network <- TF_PPI[, c(3, 4)]

  if (verbose > 0) {
    cat("\tTF network construction time:", Sys.time() - a, "\n")
  }

  # Save gene network
  store_network(network = tf_network,
                store_network = store_network,
                output_file = output_file,
                verbose = verbose)

  hummus <- add_network(hummus,
                        multiplex_name = multiplex_name,
                        network = tf_network,
                        network_name = tf_network_name,
                        verbose = verbose)

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
                                 output_file = NULL,
                                 threshold = 0.0,
                                 number_cores = 1,
                                 verbose = 1,
                                 multiplex_name = NULL,
                                 network_name = NULL) {
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
  # TODO : add other methods
  } else {
    stop(cat("Method not implemented yet, choose between GENIE3 and..",
    "that's it for now.\n but you can always compute the network",
    "independently and add it to the hummus object."))
  }

  # Save gene network
  store_network(network = gene_network,
                store_network = store_network,
                output_file = output_file,
                verbose = verbose)

  # Check if multiplex name provided
  if (is.null(multiplex_name)) {
    multiplex_name <- gene_assay
  }
  if (is.null(network_name)) {
    network_name <- paste(multiplex_name, method, sep = "_")
  }
  hummus <- add_network(hummus,
                        multiplex_name = multiplex_name,
                        network = gene_network,
                        network_name = network_name,
                        verbose = verbose)

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
    hummus,
    atac_assay = "peaks",
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    method = "cicero",
    store_network = TRUE,
    output_file = NULL,
    threshold = 0.0,
    number_cells_per_clusters = 50,
    sample_num = 100,
    seed = 2021,
    verbose = 1,
    window = 5e+05,
    reduction_method = "UMAP",
    multiplex_name = NULL,
    network_name = NULL) {

  if (method == "cicero") {
      if (!requireNamespace("cicero", quietly = TRUE)) {
        stop("Please install cicero.\n",
         "https://cole-trapnell-lab.github.io/cicero-release/docs_m3/")
      } else {

        atac_peak_network <- run_cicero_wrapper(
                                hummus,
                                atac_assay,
                                genome,
                                window,
                                number_cells_per_clusters,
                                sample_num,
                                seed,
                                verbose,
                                threshold,
                                reduction_method)
      }
  } else {
    stop(cat("Method not implemented yet, choose between Cicero and..",
    "that's it for now.\n but you can always compute the network",
    "independently and add it to the hummus object manually."))
  }
  store_network(network = atac_peak_network,
                store_network = store_network,
                output_file = output_file,
                verbose = verbose)

  if (is.null(multiplex_name)) {
    multiplex_name <- atac_assay
    }
  if (is.null(network_name)) {
    network_name <- paste0("peak_network_", method)
    }

    # Add network to hummus object
  hummus <- add_network(
    object = hummus,
    network = atac_peak_network,
    network_name = network_name,
    multiplex_name = multiplex_name,
    directed = FALSE,
    weighted = TRUE,
    verbose = verbose)

}