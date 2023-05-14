library(reticulate)
use_condaenv('base')
hummuspy <- import_from_path('hummuspy', path = '.')

g <- lapply(hummus@multilayer@multiplex, function(x) list(paste0(as.integer(x@weighted), as.integer(x@directed))))
for (multiplex in names(hummus@multilayer@multiplex)){
        names(g[[multiplex]]) =  names(hummus@multilayer@multiplex[[multiplex]]@networks)
}

b <- lapply(hummus@multilayer@bipartites, function(x) list("multiplex_right"=x@multiplex_right, "multiplex_left"=x@multiplex_left))


formatted_layers <- hummuspy$config$group_per_layer(g)

config = hummuspy$config$general_config(formatted_layers, b, bipartites_type = c('00', '00'))