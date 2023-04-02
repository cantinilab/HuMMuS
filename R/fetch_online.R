#' Title
#'
#' @param EnsDb_annotations TODO
#'
#' @return TODO
#' @export
#'
#' @examples TODO
get_genome_annotations <- function(EnsDb_annotations = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86){

  gene.range <- GetGRangesFromEnsDb(EnsDb_annotations) # Get genome annotations

  ucsc.levels <- stringr::str_replace(string=paste("chr",seqlevels(gene.range),sep=""), pattern="chrMT", replacement="chrM")
  seqlevels(gene.range) <- ucsc.levels

  return(gene.range)
}

#' Title
#'
#' @param species TODO
#'
#' @return TODO
#' @export
#'
#' @examples TODO
get_tf2motifs <- function(species = 'human'){
  #TF motifs using the union of databases: JASPAR and cis-BP included in chromVAR

  if(species=="human"){
    opts=list(collection = 'CORE',                       # Parameters for JASPAR2020
              species    = 'Homo sapiens',
              all_versions = FALSE)
    JASPAR_PWM = toPWM(getMatrixSet(JASPAR2020::JASPAR2020, opts))   # Load data from JASPAR2020
    data("human_pwms_v2")                                # Load data from chromVARmotifs
    motifs = human_pwms_v2                               # Motifs from chromVARmotifs
  }

  else if(species=="mouse"){
    opts=list(collection = 'CORE',                       # Parameters for JASPAR2020
              species    = 'Mus musculus',
              all_versions = FALSE)
    JASPAR_PWM = toPWM(getMatrixSet(JASPAR::JASPAR2020, opts))   # Load data from JASPAR2020
    data("mouse_pwms_v2")                                # Load data from chromVARmotifs
    motifs = mouse_pwms_v2                               # Motifs from chromVARmotifs
  }

  for (name in names(JASPAR_PWM)){                     # Combine motifs of JASPAR20202 and chromVARmotif
    motifs[name] = JASPAR_PWM[name]
  }
  tf2motifs <- data.frame(motif = character(),         # Initiate final TF motifs table
                          tf = character(),
                          stringsAsFactors = FALSE)
  for (i in 1:length(name(motifs))){                   # Fill TF motif table
    tfs = strsplit(name(motifs)[i], "::")[[1]]       # splitting TFs that are given as "name1::name2"
    for (tf in tfs){
      tf = strsplit(tf, "(", fixed = TRUE)[[1]][1] # only keeping <NAME> in identifier "<NAME>(var.n)"
      tf2motifs <- rbind(tf2motifs, data.frame(motif=names(motifs)[i], tf=tf))
    }
  }

  return(new("MotifsDatabase",
             tf2motifs=tf2motifs,
             motifs=motifs))
}
