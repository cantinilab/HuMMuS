#' @title Directory initialization wrapper function for the Initiate_Hummus_Object function
#' @description This function is a wrapper for the Initiate_Hummus_Object function
#' in hummus_objet.R. It creates the directories needed to store the bipartites, multiplexes, config and seed.
#' @param folder_name Directory that will store the multilayer information
  # Initialize directories
create_init_directories_wrapper <- function( folder_name ){
  multiplex_folder <- "multiplex"
  bipartite_folder <- "bipartite"
  seed_folder      <- "seed"
  config_folder    <- "config"

  dir.create(folder_name)
  dir.create(paste0(folder_name, "/", multiplex_folder))
  dir.create(paste0(folder_name, "/", bipartite_folder))
  dir.create(paste0(folder_name, "/", seed_folder))
  dir.create(paste0(folder_name, "/", config_folder))
}

#' @title Store hummus object wrapper function for the Initiate_Hummus_Object and add_network functions
#' @description This function is a wrapper for the Initiate_Hummus_Object and add_network functions
#' in hummus_objet.R, in the Initiate_Hummus_Object function it stores the initial hummus object and in the add_network one it updates the hummus object that was stored
#' @param hummus A hummus object
#' @param folder_name Directory that will store the multilayer information
  # Initialize directories
store_update_hummus_object_wrapper <- function( hummus ){
  folder_name <- hummus@multilayer_folder
  saveRDS(hummus, file.path(folder_name, "hummus_object.rds") )
}
