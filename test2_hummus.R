devtools::load_all("../hummus_package")
library(reticulate)
use_condaenv("base")
hummuspy <- import("hummuspy")


define_grn <- function(
        hummus_object,
        multiplex_names = NULL,
        bipartites_names = NULL,
        config_name = "grn_config.yml",
        config_folder = "config",
        output_f = NULL,
        tf_multiplex = "TF",
        atac_multiplex = "peaks",
        rna_multiplex = "RNA",
        multilayer_f = "multilayer",
        gene_list = NULL,
        tf_list = NULL,
        save = FALSE,
        return_df = TRUE,
        suffix_bipartites = ".tsv",
        njobs = 1) {

  print('test1')

  ##### this part should be handled with pointers
  # Check type of object
  #if (inherits(hummus_object, "multilayer")) {
   # multiplex_list <- hummus_object@multiplex
    #bipartites_list <- hummus_object@bipartites

  #} else
  if (inherits(hummus_object, "hummus_object")) {
    multiplex_list <- hummus_object@multilayer@multiplex
    bipartites_list <- hummus_object@multilayer@bipartites
  } else {
    stop("Object is not a multiplex, a multilayer nor an hummus object.")
  }

  # Check if multiplex_names is NULL
  if (is.null(multiplex_names)) {
    multiplex_names <- names(multiplex_list)
  }
  # Check if bipartites_names is NULL
  if (is.null(bipartites_names)) {
    bipartites_names <- names(bipartites_list)
  }

  # Create a named list containing the multiplexes infos
  # formatted for hummuspy config funtions
  # each element of the list is a list of the network types (weighted/directed)
  # and the name of the networks as named in the hummus object
  multiplex_dictionary <- lapply(hummus_object@multilayer@multiplex[multiplex_names],
                                function(x) list(paste0(as.integer(x@weighted), as.integer(x@directed))))
  
  print(multiplex_dictionary)
  for (multiplex in names(hummus_object@multilayer@multiplex[multiplex_names])){
    if (is.null(hummus_object@multilayer@multiplex[[multiplex]])) {
      cat('Multiplex ', multiplex, ' is NULL\n')
      next
    }

    names(multiplex_dictionary[[multiplex]]) <- names(hummus_object@multilayer@multiplex[[multiplex]]@networks)
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
  names(bipartites_dictionary) <- paste(
    names(bipartites_dictionary),
    suffix_bipartites,
    sep = "")
  
  # define GRN
  grn <- hummuspy$explore_network$define_grn_for_R(
    multilayer_f,
    multiplex_dictionary,
    bipartites_dictionary,
    gene_list = gene_list,
    tf_list = tf_list,
    config_filename = config_name,
    config_folder = config_folder,
    output_f = output_f,
    tf_multiplex = "TF",
    peak_multiplex = "peaks",
    rna_multiplex = "RNA",
    update_config = TRUE,
    save = save,
    return_df = return_df,
    njobs = njobs)

  print('test3')

  return(grn)
 }

biboo <- define_grn(
  hummus,
  multilayer_f = "b",
  njobs = 5
  )