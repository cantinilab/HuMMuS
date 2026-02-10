#' @title Directory initialization wrapper function for the Initiate_Hummus_Object function
#' @description This function is a wrapper for the Initiate_Hummus_Object function
#' in hummus_objet.R. It creates the directories needed to store the bipartites, multiplexes, config and seed.
#' @param folder_name Directory that will store the multilayer information
#' @export
#'
  # Initialize directories
create_init_directories_wrapper <- function( folder_name ){
  multiplex_folder <- "multiplex"
  bipartite_folder <- "bipartite"
  seed_folder      <- "seed"
  config_folder    <- "config"

  dir.create( folder_name, showWarnings = FALSE )
  dir.create( file.path(folder_name, multiplex_folder), showWarnings = FALSE )
  dir.create( file.path(folder_name, bipartite_folder), showWarnings = FALSE )
  dir.create( file.path(folder_name, seed_folder), showWarnings = FALSE )
  dir.create( file.path(folder_name, config_folder), showWarnings = FALSE )
}

#' @title Store hummus object wrapper function for the Initiate_Hummus_Object and add_network functions
#' @description This function is a wrapper for the Initiate_Hummus_Object and add_network functions
#' in hummus_objet.R, in the Initiate_Hummus_Object function it stores the initial hummus object and in the add_network one it updates the hummus object that was stored
#' @param hummus A hummus object
#' @export
#'
  # Initialize directories
store_update_hummus_object_wrapper <- function( hummus ){
  folder_name <- hummus@multilayer_folder
  saveRDS(hummus, file.path(folder_name, "hummus_object.rds") )
}
