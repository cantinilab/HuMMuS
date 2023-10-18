#' Format multiplex names for python hummuspy package config functions
#'
#' @param hummus_object A hummus object
#' @param multiplex_names A vector of multiplex names considered. It must be
#' a subset of the names of the multiplexes in the hummus object.
#'
#' @return A list of multiplexes names formatted for hummuspy config funtions
#' each element of the list is a list of the network types (weighted/directed)
#' and the name of the networks as named in the hummus object
#' @export
#'
#' @examples multiplexes_dictionary <- format_multiplex_names(
#'                                        hummus_object = hummus,
#'                                        multiplex_names = c("TF", "peaks"))
#'
format_multiplex_names <- function(
  hummus_object,
  multiplex_names = NULL
  ) {

    ##### this part should be handled with pointers
  # Check type of object
  # if (inherits(hummus_object, "multilayer")) {
  # multiplex_list <- hummus_object@multiplex
  # bipartites_list <- hummus_object@bipartites

  #} else
  if (inherits(hummus_object, "hummus_object")) {
    multiplex_list <- hummus_object@multilayer@multiplex
  } else {
    stop("Object is not a multilayer nor an hummus object.")
  }

  # Check if multiplex_names is NULL
  if (is.null(multiplex_names)) {
    multiplex_names <- names(multiplex_list)
  }

  # Create a named list containing the multiplexes infos
  # formatted for hummuspy config funtions
  # each element of the list is a list of the network types (weighted/directed)
  # and the name of the networks as named in the hummus object
  multiplexes_dictionary <- lapply(
    hummus_object@multilayer@multiplex[multiplex_names],
    function(x) list(paste0(as.integer(x@weighted), as.integer(x@directed))))

  # Add the names of the networks as named in the hummus object
  for (multiplex in names(hummus_object@multilayer@multiplex[multiplex_names])){
    # Check if multiplex exists in hummus object
    if (is.null(hummus_object@multilayer@multiplex[[multiplex]])) {
      cat("Multiplex ", multiplex, " is NULL\n")
      # Skip to next multiplex
      next
    }
    names(multiplexes_dictionary[[multiplex]]) <- names(
      hummus_object@multilayer@multiplex[[multiplex]]@networks)
  }
  return(multiplexes_dictionary)
  }

#' Format bipartites names for python hummuspy package config functions
#'
#' @param hummus_object A hummus object
#' @param bipartites_names A vector of bipartites names considered.
#' It must be a subset of the names of the bipartites in the hummus object.
#' @param suffix_bipartites A suffix to add to the bipartites location
#'
#' @return A list of bipartites names formatted for hummuspy config funtions
#' each element of the list is a list containing the right and left layer
#' connected by the bipartite
#' @export
#'
#' @examples bipartites_dictionary <- format_bipartites_names(
#'                                       hummus_object = hummus,
#'                                       bipartites_names = c("atac_rna",
#'                                                            "tf_peaks"))
#'
format_bipartites_names <- function(
  hummus_object,
  bipartites_names = NULL,
  suffix_bipartites = ".tsv"
  ) {

  ##### this part should be handled with pointers
  # Check type of object
  #if (inherits(hummus_object, "multilayer")) {
   # multiplex_list <- hummus_object@multiplex
    #bipartites_list <- hummus_object@bipartites

  #} else
  if (inherits(hummus_object, "hummus_object")) {
    bipartites_list <- hummus_object@multilayer@bipartites
  } else {
    stop("Object is not a multilayer nor an hummus object.")
  }

  # Check if bipartites_names is NULL
  if (is.null(bipartites_names)) {
    bipartites_names <- names(bipartites_list)
  }

  # Create a named list containing the bipartites infos
  # formatted for hummuspy config funtions
  # each element of the list is a list containing
  # the right and left layer connected by the bipartite
  bipartites_dictionary <-
              lapply(hummus_object@multilayer@bipartites[bipartites_names],
                        function(x) {
                          list("multiplex_right" = x@multiplex_right,
                               "multiplex_left" = x@multiplex_left)})
  # Add the names of the bipartites as named in the hummus object
  # and add the suffix to the names since it should indicate
  # the exact file name
  names(bipartites_dictionary) <- paste(
    names(bipartites_dictionary),
    suffix_bipartites,
    sep = "")

  # return the list
  return(bipartites_dictionary)
}

#' Define GRN from hummus object
#'
#' Calling the define_output function with output_type = 'GRN'
#'
#' @param hummus_object A hummus object
#' @param multiplex_names A vector of multiplex names considered.
#' It must be a subset of the names of the multiplexes in the hummus object.
#' @param bipartites_names A vector of bipartites names considered.
#' It must be a subset of the names of the bipartites in the hummus object.
#' @param config_name The name of the config file to be created by hummuspy
#' @param config_folder The folder where the config file will be created
#' @param tf_multiplex The name of the multiplex containing the TFs
#' @param atac_multiplex The name of the multiplex containing the ATAC-seq peaks
#' @param rna_multiplex The name of the multiplex containing the RNA-seq genes
#' @param multilayer_f The folder where the multilayer is stored
#' @param gene_list A vector of genes to be considered for the final GRN
#' (filtering is done on the genes before inferring the GRN)
#' @param tf_list A vector of TFs to be considered for the final GRN (filtering
#' is done on the TFs after inferring the GRN)
#' @param save A boolean indicating if the GRN should be saved
#' @param output_f The name of the file where the GRN should be saved
#' (if save == TRUE)
#' @param return_df A boolean indicating if the GRN should be returned as a
#' dataframe
#' @param suffix_bipartites A suffix to add to the bipartites names (to indicate
#' the exact file location)
#' @param njobs The number of jobs to be used for the computation of the GRN
#'
#' @return A dataframe containing the GRN (if return_df == TRUE)
#' @export
#'
#' @examples grn <- define_grn(hummus_object = hummus,
#'                             multilayer_f = multilayer_folder,
#'                             njobs = 5)
#'
define_grn <- function(
  hummus_object,
  multiplex_names = NULL,
  bipartites_names = NULL,
  config_name = "grn_config.yml",
  config_folder = "config",
  tf_multiplex = "TF",
  atac_multiplex = "peaks",
  rna_multiplex = "RNA",
  multilayer_f = "multilayer",
  gene_list = NULL,
  tf_list = NULL,
  save = FALSE,
  output_f = NULL,
  return_df = TRUE,
  suffix_bipartites = ".tsv",
  njobs = 1
  ) {

  grn <- define_output(
    output_type = "grn",
    hummus_object = hummus_object,
    multiplex_names = multiplex_names,
    bipartites_names = bipartites_names,
    config_name = config_name,
    config_folder = config_folder,
    tf_multiplex = tf_multiplex,
    atac_multiplex = atac_multiplex,
    rna_multiplex = rna_multiplex,
    multilayer_f = multilayer_f,
    gene_list = gene_list,
    tf_list = tf_list,
    save = save,
    output_f = output_f,
    return_df = return_df,
    suffix_bipartites = suffix_bipartites,
    njobs = njobs
  )

  # return grn
  return(grn)
 }

#' Define enhancers from hummus object
#'
#' Calling the define_output function with output_type = 'enhancers'
#'
#' @param hummus_object A hummus object
#' @param multiplex_names A vector of multiplex names considered.
#' It must be a subset of the names of the multiplexes in the hummus object.
#' @param bipartites_names A vector of bipartites names considered.
#' It must be a subset of the names of the bipartites in the hummus object.
#' @param config_name The name of the config file to be created by hummuspy
#' @param config_folder The folder where the config file will be created
#' @param tf_multiplex The name of the multiplex containing the TFs
#' @param atac_multiplex The name of the multiplex containing the ATAC-seq peaks
#' @param rna_multiplex The name of the multiplex containing the RNA-seq genes
#' @param multilayer_f The folder where the multilayer is stored
#' @param gene_list A vector of genes to be considered for the final enhancers
#' (filtering is done on the genes before inferring the enhancers)
#' @param tf_list A vector of TFs to be considered for the final enhancers
#' (filtering is done on the TFs after inferring the enhancers)
#' @param save A boolean indicating if the enhancers should be saved
#' @param output_f The name of the file where the enhancers should be saved
#' (if save == TRUE)
#' @param return_df A boolean indicating if the enhancers should be returned
#' as a dataframe
#' @param suffix_bipartites A suffix to add to the bipartites names (to indicate
#' the exact file location)
#' @param njobs The number of jobs to be used for to compute of the enhancers
#'
#' @return A dataframe containing the enhancers (if return_df == TRUE)
#' @export
#'
#' @examples enhancers <- define_enhancers(hummus_object = hummus,
#'                             multilayer_f = multilayer_folder,
#'                             njobs = 5)
#'
define_enhancers <- function(
  hummus_object,
  multiplex_names = NULL,
  bipartites_names = NULL,
  config_name = "enhancers_config.yml",
  config_folder = "config",
  tf_multiplex = "TF",
  atac_multiplex = "peaks",
  rna_multiplex = "RNA",
  multilayer_f = "multilayer",
  gene_list = NULL,
  tf_list = NULL,
  save = FALSE,
  output_f = NULL,
  return_df = TRUE,
  suffix_bipartites = ".tsv",
  njobs = 1
  ) {

  enhancers <- define_output(
    output_type = "enhancers",
    hummus_object = hummus_object,
    multiplex_names = multiplex_names,
    bipartites_names = bipartites_names,
    config_name = config_name,
    config_folder = config_folder,
    tf_multiplex = tf_multiplex,
    atac_multiplex = atac_multiplex,
    rna_multiplex = rna_multiplex,
    multilayer_f = multilayer_f,
    gene_list = gene_list,
    tf_list = tf_list,
    save = save,
    output_f = output_f,
    return_df = return_df,
    suffix_bipartites = suffix_bipartites,
    njobs = njobs
  )

  # return enhancers
  return(enhancers)
 }


#' Define binding_regions from hummus object
#'
#' Calling the define_output function with output_type = 'binding_regions'
#'
#' @param hummus_object A hummus object
#' @param multiplex_names A vector of multiplex names considered.
#' It must be a subset of the names of the multiplexes in the hummus object.
#' @param bipartites_names A vector of bipartites names considered.
#' It must be a subset of the names of the bipartites in the hummus object.
#' @param config_name The name of the config file to be created by hummuspy
#' @param config_folder The folder where the config file will be created
#' @param tf_multiplex The name of the multiplex containing the TFs
#' @param atac_multiplex The name of the multiplex containing the ATAC-seq peaks
#' @param rna_multiplex The name of the multiplex containing the RNA-seq genes
#' @param multilayer_f The folder where the multilayer is stored
#' @param gene_list A vector of genes to be considered for the final binding
#' regions (filtering is done on the genes before inferring the binding_regions)
#' @param tf_list A vector of TFs to be considered for the binding_regions
#'  (filtering is done on the TFs after inferring the binding_regions)
#' @param save A boolean indicating if the binding_regions should be saved
#' @param output_f The name of the file where the binding_regions can be saved
#' (if save == TRUE)
#' @param return_df A boolean indicating if the binding_regions should be
#' returned as a dataframe
#' @param suffix_bipartites A suffix to add to the bipartites names (to indicate
#' the exact file location)
#' @param njobs The number of jobs to be used for the computation of the binding_regions
#'
#' @return A dataframe containing the binding_regions (if return_df == TRUE)
#' @export
#'
#' @examples binding_regions <- define_binding_regions(hummus_object = hummus,
#'                             multilayer_f = multilayer_folder,
#'                             njobs = 5)
#'
define_binding_regions <- function(
  hummus_object,
  multiplex_names = NULL,
  bipartites_names = NULL,
  config_name = "binding_regions_config.yml",
  config_folder = "config",
  tf_multiplex = "TF",
  atac_multiplex = "peaks",
  rna_multiplex = "RNA",
  multilayer_f = "multilayer",
  gene_list = NULL,
  tf_list = NULL,
  save = FALSE,
  output_f = NULL,
  return_df = TRUE,
  suffix_bipartites = ".tsv",
  njobs = 1
  ) {

  binding_regions <- define_output(
    output_type = "binding_regions",
    hummus_object = hummus_object,
    multiplex_names = multiplex_names,
    bipartites_names = bipartites_names,
    config_name = config_name,
    config_folder = config_folder,
    tf_multiplex = tf_multiplex,
    atac_multiplex = atac_multiplex,
    rna_multiplex = rna_multiplex,
    multilayer_f = multilayer_f,
    gene_list = gene_list,
    tf_list = tf_list,
    save = save,
    output_f = output_f,
    return_df = return_df,
    suffix_bipartites = suffix_bipartites,
    njobs = njobs
  )

  # return binding_regions
  return(binding_regions)
 }


#' Define target genes from hummus object
#'
#' Calling the define_output function with output_type = 'target_genes'
#'
#' @param hummus_object A hummus object
#' @param multiplex_names A vector of multiplex names considered.
#' It must be a subset of the names of the multiplexes in the hummus object.
#' @param bipartites_names A vector of bipartites names considered.
#' It must be a subset of the names of the bipartites in the hummus object.
#' @param config_name The name of the config file to be created by hummuspy
#' @param config_folder The folder where the config file will be created
#' @param tf_multiplex The name of the multiplex containing the TFs
#' @param atac_multiplex The name of the multiplex containing the ATAC-seq peaks
#' @param rna_multiplex The name of the multiplex containing the RNA-seq genes
#' @param multilayer_f The folder where the multilayer is stored
#' @param gene_list A vector of genes to be considered for the target_genes
#' (filtering is done on the genes before inferring the target_genes)
#' @param tf_list A vector of TFs to be considered for the final target_genes
#' (filtering is done on the TFs after inferring the target_genes)
#' @param save A boolean indicating if the target_genes should be saved
#' @param output_f The name of the file where the target_genes should be saved
#' (if save == TRUE)
#' @param return_df A boolean indicating if the target_genes should be returned
#'  as a dataframe
#' @param suffix_bipartites A suffix to add to the bipartites names (to indicate
#' the exact file location)
#' @param njobs The number of jobs to be used to compute of the target_genes
#'
#' @return A dataframe containing the target_genes (if return_df == TRUE)
#' @export
#'
#' @examples target_genes <- define_target_genes(hummus_object = hummus,
#'                             multilayer_f = multilayer_folder,
#'                             njobs = 5)
#'
define_target_genes <- function(
  hummus_object,
  multiplex_names = NULL,
  bipartites_names = NULL,
  config_name = "target_genes_config.yml",
  config_folder = "config",
  tf_multiplex = "TF",
  atac_multiplex = "peaks",
  rna_multiplex = "RNA",
  multilayer_f = "multilayer",
  gene_list = NULL,
  tf_list = NULL,
  save = FALSE,
  output_f = NULL,
  return_df = TRUE,
  suffix_bipartites = ".tsv",
  njobs = 1
  ) {

  target_genes <- define_output(
    output_type = "target_genes",
    hummus_object = hummus_object,
    multiplex_names = multiplex_names,
    bipartites_names = bipartites_names,
    config_name = config_name,
    config_folder = config_folder,
    tf_multiplex = tf_multiplex,
    atac_multiplex = atac_multiplex,
    rna_multiplex = rna_multiplex,
    multilayer_f = multilayer_f,
    gene_list = gene_list,
    tf_list = tf_list,
    save = save,
    output_f = output_f,
    return_df = return_df,
    suffix_bipartites = suffix_bipartites,
    njobs = njobs
  )

  # return target_genes
  return(target_genes)
 }

#' @title Define output from hummus object
#'
#' @description Define output from hummus object
#'
#' @param output_type The type of output to be defined
#' @param hummus_object A hummus object
#' @param multiplex_names A vector of multiplex names considered.
#' It must be a subset of the names of the multiplexes in the hummus object.
#' @param bipartites_names A vector of bipartites names considered.
#' It must be a subset of the names of the bipartites in the hummus object.
#' @param config_name The name of the config file to be created by hummuspy
#' @param config_folder The folder where the config file will be created
#' @param tf_multiplex The name of the multiplex containing the TFs
#' @param atac_multiplex The name of the multiplex containing the ATAC-seq peaks
#' @param rna_multiplex The name of the multiplex containing the RNA-seq genes
#' @param multilayer_f The folder where the multilayer is stored
#' @param gene_list A vector of genes to be considered for the target_genes
#' (filtering is done on the genes before inferring the target_genes)
#' @param tf_list A vector of TFs to be considered for the final target_genes
#' (filtering is done on the TFs after inferring the target_genes)
#' @param save A boolean indicating if the target_genes should be saved
#' @param output_f The name of the file where the target_genes should be saved
#' (if save == TRUE)
#' @param return_df A boolean indicating if the target_genes should be returned
#'  as a dataframe
#' @param suffix_bipartites A suffix to add to the bipartites names (to indicate
#' the exact file location)
#' @param njobs The number of jobs to be used to compute of the target_genes
#'
#' @return A dataframe containing the target_genes (if return_df == TRUE)
#' @export
#'
#' @examples target_genes <- define_output('grn', hummus_object = hummus)
define_output <- function(
  output_type,
  hummus_object,
  multiplex_names = NULL,
  bipartites_names = NULL,
  config_name = "config.yml",
  config_folder = "config",
  tf_multiplex = "TF",
  atac_multiplex = "peaks",
  rna_multiplex = "RNA",
  multilayer_f = "multilayer",
  gene_list = NULL,
  tf_list = NULL,
  save = FALSE,
  output_f = NULL,
  return_df = TRUE,
  suffix_bipartites = ".tsv",
  njobs = 1
  ) {

  # Check if hummuspy is installed and import it
  hummuspy <- tryCatch({
    reticulate::import("hummuspy")
    }, error = function(err) {
      stop("hummuspy package not found. Make sure that Reticulate \
      is pointing to the right Python binary.")
      }
  )
  # Format multiplexes names
  multiplexes_dictionary <- format_multiplex_names(
    hummus_object,
    multiplex_names = multiplex_names)
  # Format bipartites names
  bipartites_dictionary <- format_bipartites_names(
    hummus_object,
    bipartites_names = bipartites_names,
    suffix_bipartites = suffix_bipartites)

  # define target_genes with hummuspy function
  output <- hummuspy$explore_network$get_output_from_dicts(
    output_request = output_type,
    multilayer_f = multilayer_f,
    multiplexes_list = multiplexes_dictionary,
    bipartites_list = bipartites_dictionary,
    gene_list = gene_list,
    tf_list = tf_list,
    config_filename = config_name,
    config_folder = config_folder,
    output_f = output_f,
    tf_multiplex = tf_multiplex,
    peak_multiplex = atac_multiplex,
    rna_multiplex = rna_multiplex,
    update_config = TRUE,
    save = save,
    return_df = return_df,
    njobs = njobs)

  # return target_genes
  return(output)
 }

#' @title Define general config file for hummuspy
#' 
#' @description Define general config file for hummuspy
#' 
#' @param hummus_object A hummus object
#' @param multiplex_names A vector of multiplex names considered.
#'  It must be a subset of the names of the multiplexes in the hummus object, or NULL
#'  if all multiplexes should be considered.
#' @param bipartites_names A vector of bipartites names considered.
#' It must be a subset of the names of the bipartites in the hummus object, or NULL
#' if all bipartites should be considered.
#' @param folder_multiplexes The folder where the multiplexes are stored
#' @param folder_bipartites The folder where the bipartites are stored
#' @param seed_path The path to the seed file
#' @param suffix_bipartites A suffix to add to the bipartites names (to indicate
#' the exact file name)
#' @param self_loops A boolean indicating if self loops should be considered.
#' @param restart_proba The restart probability for the random walk (default = 0.7)
#' @param save_configfile A boolean indicating if the config file should be saved
#' @param config_name The name of the config file to be created by hummuspy
#' @param config_folder The folder where the config file will be created (inside multilayer_f)
#' @param multilayer_f The folder where the multilayer is stored
#' 
#' @return A config file for hummuspy
#' @export
#' 
define_general_config <- function(
  hummus_object,
  multiplex_names = NULL,
  bipartites_names = NULL,
  folder_multiplexes = "multiplex",
  folder_bipartites = "bipartites",
  seed_path = 'seed/seeds.txt',
  suffix = ".tsv",
  self_loops = FALSE,
  restart_proba = 0.7,
  save_configfile = FALSE,
  config_name = "config.yml",
  config_folder = "config",
  multilayer_f = "multilayer"
  ) {

  # Check if hummuspy is installed and import it
  hummuspy <- tryCatch({
    reticulate::import("hummuspy")
    }, error = function(err) {
      stop("hummuspy package not found. Make sure that Reticulate \
      is pointing to the right Python binary.")
      }
  )
  # Format multiplexes names
  multiplexes_dictionary <- format_multiplex_names(
    hummus_object,
    multiplex_names = multiplex_names)
  # Format bipartites names
  bipartites_dictionary <- format_bipartites_names(
    hummus_object,
    bipartites_names = bipartites_names,
    suffix_bipartites = suffix_bipartites)

  self_loops <- as.integer(self_loops)

  if (save_configfile == TRUE) {
    config_filename <- file.path(mutlilayer_f, config_folder, config_name)
  } else {
    config_filename <- NULL
  }

  # define target_genes with hummuspy function
  config <- hummuspy$config$general_config(
    multiplexes = multiplexes_dictionary,
    bipartites = bipartites_dictionary,
    folder_multiplexes = folder_multiplexes,
    folder_bipartites = folder_bipartites,
    seed_path = seed_path,
    self_loops = self_loops,
    restart_prob = restart_proba,
    config_filename = config_filename,
    save_configfile = save_configfile,
    suffix = suffix)

  return(config)
 }
