#' @title Motifs database class
#' MotifsDatabase object stores motifs(PFM matrices)
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
                             tfs = "NULL"
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
      paste("Multiplex of ", length(object@networks),
      " networks with", length(object@features), "features.\n",
      "Networks names: ", paste(names(object@networks), collapse = ", "))
      )
  })

add_network <- function(
  object,
  network,
  network_name,
  multiplex_name = NULL,
  directed = FALSE,
  weighted = FALSE) {

  #Check type of object
  if (inherits(object, "multiplex")) {
    multiplex <- object
  } else if (inherits(object, "multilayer")) {
    if (is.null(multiplex_name)) {
      stop("You need to specify the multiplex name.")
    }
    multiplex <- object@multiplex[[multiplex_name]]
  } else if (inherits(object, "hummus_object")) {
    if (is.null(multiplex_name)) {
      stop("You need to specify the multiplex name.")
    }
    multiplex <- object@multilayer@multiplex[[multiplex_name]]
  } else {
    stop("Object is not a multiplex, a multilayer nor an hummus object.")
  }
  # Check if network name already exists in the multiplex
  if (network_name %in% names(multiplex@networks)) {
    stop("Network name already exists")
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
    object@multiplex <- multiplex
    return(object)
  } else if (inherits(object, "hummus_object")) {
    object@multilayer@multiplex <- multiplex
    return(object)
  }
}

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

multilayer <- setClass(Class = "multilayer",
                       slots = c(
                        "bipartites" = "list", # Bipartite networks
                        "multiplex" = "list", # Multiplex networks
                        "config" = "list" # Parameters to create the hmln
                        )                 # representation of a yaml file
                      )
setMethod("show", "multilayer",
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
#' The SeuratPlus object is an extended \code{Seurat} object
#' for the storage and analysis of a heterogeneous multilayer network
#'
#' @slot multilayer
#'
#' @name hummus_object-class
#' @rdname hummus_object-class
#' @exportClass hummus_object
#' @concept assay
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

save_multilayer <- function(
    hummus,
    folder_name,
    verbose = TRUE,
    suffix = ".tsv") {

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