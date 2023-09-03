#' @importFrom methods setClass
#' @importClassesFrom Signac Motif
#' @importClassesFrom SeuratObject Seurat
#' @importClassesFrom TFBSTools PWMatrixList
NULL


#' @title Motifs database class
#'
#' @description MotifsDatabase object stores motifs(PFM matrices)
#' and tf2motifs (TF to motifs names mapping) data.
#'
#' @slot motifs (TFBSTools::PWMatrixList) - PFM matrices.
#' @slot tf2motifs (data.frame) - TF to motif names mapping. Columns: motif, tf.
#'
#' @name  motifs_db-class
#' @rdname motifs_db-class
#' @exportClass motifs_db
motifs_db <- setClass("motifs_db",
                           representation(
                             motifs = "PWMatrixList",
                             tf2motifs = "data.frame",
                             tfs = "character"
                           ))
setMethod("show", "motifs_db",
  function(object) {
    cat(
      paste("Motifs database object with :\n- ",
          length(object@motifs), "motifs\n- ",
          length(unique(object@tf2motifs$tf)), " TFs\n- ",
          nrow(object@tf2motifs), "TF to motif names mapping"
          )
      )
  })


#' @title Multiplex class
#' @description Multiplex object stores a list of networks, a list of features and
#' a list of logicals indicating if the network is directed or weighted.
#' @slot networks (list) - List of networks.
#' @slot features (vector) - Vector of features.
#' @slot directed (list) - List of logical indicating if networks are directed.
#' @slot weighted (list) - List of logical indicating if networks are weighted.
#' 
#' @name multiplex-class
#' @rdname multiplex-class
#' @exportClass multiplex
multiplex <- setClass(Class = "multiplex",
                       slots = c(
                         "networks" = "list", # List of networks
                         "features" = "vector", # Vector of features
                         "directed" = "list", # Logical indicating
                                      # if the network is directed
                         "weighted" = "list" # Logical indicating
                                      # if the network is weighted
                         # "network_names" = "vector" # Vector of network names
                        )
                      )

setMethod("show", "multiplex",
  function(object) {
    cat(
      # Reprensentation of the multiplex object
      # with the number of networks and features, and the list of network names
      paste("Multiplex of ", length(object@networks),
      " networks with", length(object@features), "features.\n",
      "Networks names: ", paste(names(object@networks), collapse = ", "))
      )
  })


#' @title Bipartite class
#'
#' @description Bipartite object stores a bipartite network (edge list) and the names of the
#'  left and right features' multiplexes.
#' @slot network (data.frame) - Bipartite network (edge list)
#' @slot multiplex_left (character) - Left features' multiplex
#' @slot multiplex_right (character) - Right features' multiplex
#'
#' @name bipartite-class
#' @rdname bipartite-class
#' @exportClass bipartite
#'
#' @examples bipartite <- bipartite(
#'                           network = bipartite_network,
#'                          multiplex_left = "RNA",
#'                         multiplex_right = "peaks")
#'
bipartite <- setClass(Class = "bipartite",
                       slots = c(
                      "network" = "data.frame", # Bipartite network (edge list)
                      "multiplex_left" = "character", # left features' multiplex
                      "multiplex_right" = "character" # right features multiplex
                        )
                      )

setMethod("show", "bipartite",
  function(object) {
    cat(
      paste("Bipartite network of ", nrow(object@network), " edges.\n",
      "Multiplexes names: ", object@multiplex_left,
      " and ", object@multiplex_right, "\n")
      )
  })

#' @title Multilayer class
#'
#' @description Multilayer object stores a list of bipartite networks and a
#'  list of multiplex networks. It can also stores a config list to create a
#'  yaml file, which is used to parametrize the random walk with restart to
#' explore the multilayer.
#'
#' @slot bipartites (list) - List of bipartite networks
#' @slot multiplex (list) - List of multiplex networks
#' @slot config (list) - List of parameters to parametrize the random walk with
#' restart to explore the multilayer
#'
#' @name multilayer-class
#' @rdname multilayer-class
#' @exportClass multilayer
#'
multilayer <- setClass(Class = "multilayer",
                       slots = c(
                        "bipartites" = "list", # Bipartite networks
                        "multiplex" = "list", # Multiplex networks
                        "config" = "list" # Parameters to create the hmln
                        )                 # representation of a yaml file
                      )

setMethod("show", "multilayer",
  # Representation of the multilayer object with the number of bipartite and
  # multiplex networks, and the list of bipartite names and multiplex names
  function(object) {
    cat(
      paste("Multilayer network containing ",
      length(object@bipartites), " bipartite networks and ",
      length(object@multiplex), " multiplex networks.\n",
      "\n- Multiplex names: ", paste(names(object@multiplex),
                                          collapse = ", "),
      "\n- Bipartite names: ", paste(names(object@bipartites),
                                          collapse = ", "), "\n"
      )
    )
  })


#' The hummus_object class
#'
#' The hummus_object object is an extended \code{Seurat} object
#' for the storage and analysis of a heterogeneous multilayer network
#'
#' @slot multilayer (multilayer) - Multilayer object
#' @slot motifs_db (motifs_db) - Motifs database
#' @slot assay (list) - List of assays
#'
#' @name hummus_object-class
#' @rdname hummus_object-class
#' @exportClass hummus_object
#' @export
#'
#' @examples hummus_object <- hummus_object(seurat_object)
#'
hummus_object <- setClass(
    Class = "hummus_object",
    contains = "Seurat",
    slots = list(
        "multilayer" = "multilayer",
        "motifs_db" = "motifs_db"
    )
)

setMethod("show", "hummus_object",
  function(object) {
    object <- SeuratObject::UpdateSlots(object = object)
    assays <- SeuratObject::FilterObjects(object = object,
                                          classes.keep = "Assay")
    nfeatures <- sum(vapply(
      X = assays,
      FUN = function(x) {
        return(nrow(x = object[[x]]))
      },
      FUN.VALUE = integer(length = 1L)
    ))
    num.assays <- length(x = assays)

    cat("Hummus object containing a multilayer object :\n")
    show(object@multilayer)
    cat('\n\nAnd a Seurat object :\n\n')
    cat(
      nfeatures,
      "features across",
      ncol(x = object),
      "samples within",
      num.assays,
      ifelse(test = num.assays == 1, yes = "assay", no = "assays"),
      "\n"
    )
    cat(
      "Active assay:",
      SeuratObject::DefaultAssay(object = object),
      paste0('(', nrow(x = object), " features, ",
      length(x = SeuratObject::VariableFeatures(object = object))," variable features)")
    )
    other.assays <- assays[assays != SeuratObject::DefaultAssay(object = object)]
    if (length(x = other.assays) > 0) {
      cat(
        '\n',
        length(x = other.assays),
        'other',
        ifelse(test = length(x = other.assays) == 1, yes = 'assay', no = 'assays'),
        'present:',
        strwrap(x = paste(other.assays, collapse = ', '))
      )
    }
    reductions <- SeuratObject::FilterObjects(object = object, classes.keep = 'DimReduc')
    if (length(x = reductions) > 0) {
      cat(
        '\n',
        length(x = reductions),
        'dimensional',
        ifelse(test = length(x = reductions) == 1, yes = 'reduction', no = 'reductions'),
        'calculated:',
        strwrap(x = paste(reductions, collapse = ', '))
      )
    }
    fovs <- SeuratObject::FilterObjects(object = object, classes.keep = 'FOV')
    if (length(x = fovs)) {
      cat(
        '\n',
        length(x = fovs),
        'spatial',
        ifelse(test = length(x = fovs) == 1L, yes = 'field', no = 'fields'),
        'of view present:',
        strwrap(x = paste(fovs, sep = ', '))
      )
    }
    images <- SeuratObject::FilterObjects(object = object, classes.keep = 'SpatialImage')
    images <- setdiff(x = images, y = fovs)
    if (length(x = images)) {
      cat(
        '\n',
        length(x = images),
        ifelse(test = length(x = images) == 1L, yes = 'image', no = 'images'),
        'present:',
        strwrap(x = paste(images, collapse = ', '))
      )
    }
    cat('\n')
  }
)


#' @title Save multilayer object files in a hierarchical structure on disk
#'
#' @description Save multilayer files from a hummus_object
#' in a hierarchical structure on disk, inside a folder specified through
#'  folder_name
#'
#' @param hummus A hummus object
#' @param folder_name The name of the folder to save the multilayer
#' @param verbose (integer) - Display function messages. Set to 0 for no
#'  message displayed, >= 1 for more details.
#' @param suffix The suffix of the files to save. Default: ".tsv"
#'
#' @return Nothing, but create a folder containing the multilayer object files
#' @export
#'
#' @examples folder_name = "multilayer"
#' save_multilayer(hummus = hummus, folder_name = "multilayer")
#'
save_multilayer <- function(
    hummus,
    folder_name,
    verbose = TRUE,
    suffix = ".tsv"
    ) {

  multiplex_folder <- "multiplex"
  bipartite_folder <- "bipartite"
  seed_folder      <- "seed"
  config_folder    <- "config"

  dir.create(folder_name)
  dir.create(paste0(folder_name, "/", multiplex_folder))
  dir.create(paste0(folder_name, "/", bipartite_folder))
  dir.create(paste0(folder_name, "/", seed_folder))
  dir.create(paste0(folder_name, "/", config_folder))

  # For each multiplex, create a subfolder of multiplex, 
  # and save its networks inside
  for (multiplex_name in names(hummus@multilayer@multiplex)){
    dir.create(paste0(folder_name, "/", multiplex_folder, "/", multiplex_name))
    print(hummus@multilayer@multiplex[[multiplex_name]])
    for (network_name in names(hummus@multilayer@multiplex[[multiplex_name]]@networks)){
      print(paste(multiplex_name, network_name))
      write.table(hummus@multilayer@multiplex[[multiplex_name]]@networks[[network_name]],
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t",
             file = paste0(folder_name, "/",
                           multiplex_folder, "/",
                           multiplex_name, "/", network_name, suffix))
    }
  }
  # save bipartite networks
  for (bipartite in names(hummus@multilayer@bipartites)){
      write.table(hummus@multilayer@bipartites[[bipartite]]@network, sep = "\t",
                  col.names = FALSE, row.names = FALSE, quote = FALSE,
                  file = paste0(folder_name, "/",
                               bipartite_folder, "/",
                               bipartite, ".tsv"))
  }
}


#' @title Add a network to a multiplex, a multilayer or an hummus object
#'
#' @description Add a network to a multiplex, a multilayer or an hummus object
#'
#' @param object A multiplex, a multilayer or an hummus object
#' @param network A network (edge list)
#' @param network_name The name of the network
#' @param multiplex_name The name of the multiplex. Default: NULL if object is a
#' multiplex already only
#' @param directed Logical indicating if the network is directed. Default: FALSE
#' @param weighted Logical indicating if the network is weighted. Default: FALSE
#' @param verbose (integer) - Display function messages. Set to 0 for no
#' message displayed, >= 1 for more details.
#'
#' @return A multiplex, a multilayer or an hummus object with the added network
#' @export
#'
#' @examples hummus <- add_network(
#'                            object = hummus,
#'                            network = atac_peak_network,
#'                            network_name = network_name,
#'                            multiplex_name = multiplex_name,
#'                            weighted = TRUE,
#'                            directed = FALSE)
#' 
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
    if (!(multiplex_name %in% names(object@multiplex))) {
      if (verbose > 0) {
        cat("\tCreating new multiplex : ", multiplex_name, "\n")
      }
      # Create new multiplex if not
      object@multiplex[[multiplex_name]] <- new("multiplex")
    }
    # Get working multiplex
    multiplex <- object@multiplex[[multiplex_name]]
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


#' @title Wrapper function to save a network or not
#'
#' @description Wrapper function to save a network or not in a file according
#' to the store_network parameter. If store_network is TRUE, the network is
#' saved in the output_file.
#'
#' @param network A network (edge list)
#' @param store_network Logical indicating if the network should be saved
#' @param output_file The name of the file to save the network
#' @param verbose (integer) - Display function messages. Set to 0 for no
#' message displayed, >= 1 for more details.
#'
#' @return Nothing, but save the network in a file if store_network is TRUE
#' @export
#'
#' @examples network <- read.table("network.tsv", header = TRUE, sep = "\t")
#'           store_network(network = network,
#'               store_network = TRUE,
#'               output_file = "network.tsv",
#'               verbose = 1)
#'
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