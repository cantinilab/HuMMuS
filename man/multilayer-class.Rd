% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hummus_objet.R
\docType{class}
\name{multilayer-class}
\alias{multilayer-class}
\alias{multilayer}
\title{Multilayer class}
\description{
Multilayer object stores a list of bipartite networks and a
list of multiplex networks. It can also stores a config list to create a
yaml file, which is used to parametrize the random walk with restart to
explore the multilayer.
}
\section{Slots}{

\describe{
\item{\code{bipartites}}{(list) - List of bipartite networks}

\item{\code{multiplex}}{(list) - List of multiplex networks}

\item{\code{config}}{(list) - List of parameters to parametrize the random walk with
restart to explore the multilayer}
}}

