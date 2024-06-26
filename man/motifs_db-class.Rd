% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hummus_objet.R
\docType{class}
\name{motifs_db-class}
\alias{motifs_db-class}
\alias{motifs_db}
\title{Motifs database class}
\description{
MotifsDatabase object stores motifs(PFM matrices)
and tf2motifs (TF to motifs names mapping) data.
}
\section{Slots}{

\describe{
\item{\code{motifs}}{(TFBSTools::PWMatrixList) - PFM matrices.}

\item{\code{tf2motifs}}{(data.frame) - TF to motif names mapping. Columns: motif, tf.}
}}

