#' Extract multilayer and bipartites from hummus object
#'
#' @description Extract multilayer and bipartites from hummus object,
#'  in order to use them in the python script building the config file.
#'
#' param hummus_object: hummus object
#' param layers: list of layers to extract if not all
#' param bipartites: list of bipartites to extract if not all
#'
#'
extract_names <- function(
    hummus_object,
    layers = NULL,
    bipartites = NULL) {

    # Extract layers
    if (is.null(layers)) {
        layers <- names(hummus_object@multilayer@multiplex)
    }
    if (is.null(bipartites)) {
        bipartites <- names(hummus_object@multilayer@bipartites)
    }

}