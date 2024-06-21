#' Compute TF network and add it to hummus object
#'
#' Compute a protein-protein interaction layer from Omnipath request that will represent tf cooperativity.
#' This network is the top-layer of HuMMuS multilayer.
#'
#' @param hummus (Hummus_Object) - Hummus object
#' @param organism (integer)  - Specie identifier from Omnipath to fetch
#' specific interactions
#' @param tfs vector(character) - List of tfs consider. If NA, tfs are extracted
#' from the hummus object with get_tfs function.
#' @param gene_assay (character) - Name of the assay to get tfs from if tfs is
#' not provided. If NULL, all TFs with motifs in the hummus object are used.
#' @param method (character) - Method used to infer network edges.
#' * \code{'Omnipath'} - Use Omnipath to infer tf-tf networks.
#' * \code{'NULL'} - A fake connected network is computed.
#' * \code{'Other method'} - TO DO.
#' @param store_network (bool) - Save the network directly (\code{TRUE},
#'  default) or return without saving on disk (\code{FALSE}).
#' @param output_file (character) - Name of the output_file
#' (if store_network == \code{TRUE}).
#' @param source_target ('AND'|'OR') - Fetch only the interactions involving
#' two considered tfs (\code{'AND', default}), or one considered tfs and any 
#' other element (\code{'OR'})
#' @param multiplex_name (character) - Name of the multiplex to add the network
#'  to. Default is \code{'TF'}.
#' @param tf_network_name (character) - Name of the network in the multiplex to
#' add the network to. Default is \code{'TF_network'}.
#' @param verbose (integer) - Display function messages. Set to 0 for no message
#'  displayed, >= 1 for more details.
#'
#' @return (Hummus_Object) - Return hummus object with the new network added.
#' @export
#'
#' @examples hummus <- compute_tf_network(hummus,
#'                                        gene_assay = "RNA",
#'                                        verbose = 1)
compute_tf_network <- function(
  hummus, # Hummus object
  organism = 9606, # Human by default
  tfs = NA, # List of tfs considered.
  gene_assay = NULL, # Name of the assay to get tfs from
                     # if tfs is not provided
  method = NULL, # Method used to infer network edges.
                      # * 'Omnipath' - Use Omnipath to infer tf-tf networks.
                      # * 'NULL' - A fake connected network is computed.
                      # * 'Other method' - TO DO.
  store_network = FALSE, # Save the network on disk (TRUE, default)
  output_file = NULL, # Name of the output_file (if store_network == TRUE)
  source_target = "AND", # 'AND' | 'OR'
  multiplex_name = "TF", # Name of the multiplex to add the network to
  tf_network_name = "TF_network", # Name of the network in the multiplex
  verbose = 1
  ) {

  a <- Sys.time()
  # Check if method is implemented
  if (is.null(method)) {
        tf_network <- run_tf_null_wrapper(
          hummus = hummus,
          organism = organism,
          tfs = tfs,
          gene_assay = gene_assay,
          verbose)
  } else if (method == "Omnipath") {
      if (!requireNamespace("OmnipathR", quietly = TRUE)) {
        stop("Please install Omnipath.\n",
          "github.com/saezlab/OmnipathR")
      } else {
        # infer network with cicero
        tf_network <- run_omnipath_wrapper(
          hummus = hummus,
          organism = organism,
          tfs = tfs,
          gene_assay = gene_assay,
          source_target = source_target,
          verbose = verbose)
      }
  }  else {
    stop(cat("Method not implemented yet, choose between Omnipath and NULL..",
    "that's it for now.\n But you can always compute the network",
    "independently and add it to the hummus object manually !"))
  }
  if (verbose > 0) {
    cat("TF network construction time:", Sys.time() - a)
  }

  # Save gene network
  store_network(network = tf_network,
                store_network = store_network,
                output_file = output_file,
                verbose = verbose)

  # Add network to hummus object
  hummus <- add_network(hummus,
                        multiplex_name = multiplex_name,
                        network = tf_network,
                        network_name = tf_network_name,
                        weighted = FALSE, # PPI could be weighted,
                                          # could be added later
                        directed = FALSE, # PPI are not directed
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
#' @param hummus (Hummus_Object) - Hummus object
#' @param gene_assay (character) - Name of the assay containing the gene
#'  expression data.
#' @param tfs vector(character) - List of tfs considered. If NULL, all TFs with
#' motifs in the hummus object are used.
#' @param method (character) - Method used to infer network edges.
#' * \code{'Genie3'} - Use tree random forest to infer regulatory networks.
#' * \code{'Other method'} - TO DO.
#' @param multiplex_name (character) - Name of the multiplex to add the network
#' to. Default is \code{'RNA'}.
#' @param network_name (character) - Name of the network in the multiplex to
#' add the network to. Default is \code{'RNA_network'}.
#' @param store_network (bool) - Save the network directly (\code{TRUE},
#'  default) or return without saving on disk (\code{FALSE}).
#' @param output_file (character) - Name of the output_file
#'  (if store_network == \code{TRUE}).
#' @param threshold (interger, default 0) - Minimal threshold
#'  to select tf-gene edges.
#' @param number_cores (interger, default 1) - Number of thread that should be
#'  used for the parallelizable methods.
#' @param verbose (integer) - Display function messages. Set to 0 for no
#'  message displayed, >= 1 for more details.
#'
#' @return (data.frame) - Return list of network interactions between genes
#' @export
#'
#' @examples hummus <- compute_gene_network(
#'                                hummus,
#'                                gene_assay = "RNA",
#'                                method = "GENIE3",
#'                                verbose = 1,
#'                                number_cores = 8,
#'                                store_network = FALSE)
#' 
compute_gene_network <- function(
  hummus,
  gene_assay = "RNA",
  tfs = NULL,
  method = "GENIE3",
  multiplex_name = NULL,
  network_name = NULL,
  store_network = FALSE,
  output_file = NULL,
  threshold = 0.0,
  number_cores = 1,
  verbose = 1
  ) {

  a <- Sys.time()
  # Check if method is implemented
  if (method == "GENIE3") {
    if (verbose > 0) {
      cat("Computing gene network with ", method, " ...\n")
    }
    # Get tfs list
    if (verbose > 0 && is.null(tfs)) {
      cat("\tNo TFs list provided, fetching from hummus object...\n")
      tfs <- get_tfs(hummus = hummus,
            assay = gene_assay,
            store_tfs = FALSE,
            output_file = NULL,
            verbose = verbose)
    }

    # infer network
    weightMat <- GENIE3::GENIE3(as.matrix(hummus@assays[[gene_assay]]@counts),
                               regulators = tfs,
                               nCores = number_cores)
    # get edge list
    linkList <- GENIE3::getLinkList(weightMat)
    gene_network <- linkList[which(linkList$weight > threshold), ]
 
  # TODO : add other methods
  } else {
    stop(cat("Method not implemented yet, choose between GENIE3 and..",
    "that's it for now.\n but you can always compute the network",
    "independently and add it to the hummus object."))
  }
  if (verbose > 0) {
      cat("\tGene network construction time:", Sys.time() - a, "\n")
  }

  # Save gene network
  store_network(network = gene_network,
                store_network = store_network,
                output_file = output_file,
                verbose = verbose)

  # If no multiplex name provided, use assay name
  if (is.null(multiplex_name)) {
    multiplex_name <- gene_assay
  }
  # If no network name provided, use method name + assay name
  if (is.null(network_name)) {
    network_name <- paste(multiplex_name, method, sep = "_")
  }
  # Add network to hummus object
  hummus <- add_network(hummus,
                        multiplex_name = multiplex_name,
                        network = gene_network,
                        network_name = network_name,
                        weighted = TRUE,
                        directed = FALSE,
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
#' @param hummus (Hummus_Object) - Hummus object
#' @param atac_assay (character) - Name of the assay containing the atac
#'  peaks data.
#' @param genome (BSgenome) - Genome used to compute the distance between peaks.
#' @param method (character) - Method used to infer network edges.
#' * \code{'cicero'} - Use cicero to infer regulatory networks.
#' * \code{'Other method'} - TO DO.
#' @param multiplex_name (character) - Name of the multiplex to add the network
#' to. Default is \code{'peaks'}.
#' @param network_name (character) - Name of the network in the multiplex to
#' add the network to. Default is \code{'peak_network'}.
#' @param store_network (bool) - Save the network directly (\code{TRUE},
#'  default) or return without saving on disk (\code{FALSE}).
#' @param output_file (character) - Name of the output_file
#'  (if store_network == \code{TRUE}).
#' @param threshold (interger, default 0) - Minimal threshold to select tf-gene
#'  edges.
#' @param number_cells_per_clusters (integer) - Number of cells grouped by
#'  territory to define pseudocells
#' @param sample_num (integer | Cicero) - Number of pseudocells to sample from
#' each territory. Default is 100.
#' @param seed (integer | Cicero) - Seed used to sample pseudocells. Default is
#' 2025
#' @param verbose (integer) - Display function messages. Set to 0 for no
#'  message displayed, >= 1 for more details.
#' @param window (integer) - Size of window to consider potential
#'  cis-regulatory cooperations between peaks. Default is 500K base pairs.
#' @param reduction_method (character | Cicero) - Method used to reduce dimensionality
#' of the data to identify territories. Default is \code{'UMAP'}.
#'
#' @return (data.frame) - Return list of network interactions between peaks
#' @export
#'
#' @examples hummus <- compute_atac_peak_network(hummus)
#' 
compute_atac_peak_network <- function(
    hummus,
    atac_assay = "peaks",
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    method = "cicero",
    multiplex_name = NULL,
    network_name = NULL,
    store_network = FALSE,
    output_file = NULL,
    threshold = 0.0,
    number_cells_per_clusters = 50,
    sample_num = 100,
    seed = 2025,
    verbose = 1,
    window = 5e+05,
    reduction_method = "UMAP") {
    
  a <- Sys.time()
  # Check if method is implemented
  if (method == "cicero") {
      if (!requireNamespace("cicero", quietly = TRUE)) {
        stop("Please install cicero.\n",
         "https://cole-trapnell-lab.github.io/cicero-release/docs_m3/")
      } else {
        # infer network with cicero
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
  if (verbose > 0) {
    cat("Peak network construction time:", Sys.time() - a)
  }
  # Save peak network
  store_network(network = atac_peak_network,
                store_network = store_network,
                output_file = output_file,
                verbose = verbose)
  # If no multiplex name provided, use assay name
  if (is.null(multiplex_name)) {
    multiplex_name <- atac_assay
    }
  # If no network name provided, use method name + assay name
  if (is.null(network_name)) {
    network_name <- paste0("peak_network_", method)
    }

  # Add network to hummus object
  hummus <- add_network(
    object = hummus,
    network = atac_peak_network,
    network_name = network_name,
    multiplex_name = multiplex_name,
    weighted = TRUE,
    directed = FALSE,
    verbose = verbose)

}
