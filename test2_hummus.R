library(reticulate)
use_condaenv('base')
hummuspy <- import_from_path('hummuspy', path = '.')

eta <-list(0,1,0)
lamb <- list(list(1/3,1/3,1/3),
             list("1/3","1/3","1/3"),
             list("1/3","1/3","1/3"))

create_config <- function(
        hummus_object,
        multiplex_names = NULL,
        bipartites_names = NULL,
        eta = NULL,
        lamb = NULL) {

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
  bipartites_dictionary <- lapply(hummus_object@multilayer@bipartites[bipartites_names],
                                 function(x) list("multiplex_right"=x@multiplex_right, "multiplex_left"=x@multiplex_left))


  #  Reformatting the list to be used in hummuspy config functions
  formatted_layers <- hummuspy$config$group_per_layer(multiplex_dictionary)
  print(formatted_layers)

  # Initialize the config object with hummus object
  # multiple layers and bipartites infos
  config <- hummuspy$config$general_config(formatted_layers, bipartites_dictionary, bipartites_type = c('00', '00'))

  config <- hummuspy$config$setup_proba_config(config, eta, lamb)
  # CREATE A SAVE CONDITION
  hummuspy$config$save_config(config, "a/config.yaml")
}

create_config(
        hummus,
        eta = NULL,
        lamb = NULL)