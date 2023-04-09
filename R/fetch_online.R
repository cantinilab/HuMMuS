#'  Fetch online genome annotations from Ensembldb database
#'
#' @param EnsDb_annotations (EndsDb object) - Ensembldb database (default: EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
#'
#' @return gene_range (GRanges object) - Genome annotations
#' @export
#'
#' @examples gene_range = get_genome_annotations(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
get_genome_annotations <- function(
  ensdb_annotations = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86) {

  gene_range <- Signac::GetGRangesFromEnsDb(ensdb_annotations) # Get genome annotations

  ucsc.levels <- stringr::str_replace(
    string = paste("chr", seqlevels(gene_range), sep = ""),
    pattern = "chrMT",
    replacement = "chrM")

  seqlevels(gene_range) <- ucsc.levels

  return(gene_range)
}

#' Fetch online TF motifs from JASPAR2020 and chromVARmotifs
#'
#' @param species TODO
#'
#' @return TODO
#' @export
#'
#' @examples TODO
get_tf2motifs <- function(species = "human") {
  #TF motifs using the union of databases: JASPAR and cis-BP
  # included in chromVAR
  getMatrixSet <- TFBSTools::getMatrixSet
  if (species == "human") {
    # Parameters for JASPAR2020
    opts=list(collection = "CORE",
              species    = "Homo sapiens",
              all_versions = FALSE)
    JASPAR_PWM <- TFBSTools::toPWM(getMatrixSet(JASPAR2020::JASPAR2020, opts))
    # Load data from JASPAR2020
    # Load data from chromVARmotifs
    # Original data accessible at https://github.com/GreenleafLab/chromVARmotifs
    data("human_pwms_v2")
    # Load data from chromVARmotifs
    motifs <- human_pwms_v2
    # Motifs from chromVARmotifs
  } else if (species == "mouse") {
    # Parameters for JASPAR2020
    opts=list(collection = "CORE",
              species    = "Mus musculus",
              all_versions = FALSE)
    JASPAR_PWM <- TFBSTools::toPWM(getMatrixSet(JASPAR::JASPAR2020, opts))
    # Load data from JASPAR2020
    data("mouse_pwms_v2")
    # Load data from chromVARmotifs
    # Original data accessible at https://github.com/GreenleafLab/chromVARmotifs
    motifs <- mouse_pwms_v2
    # Motifs from chromVARmotifs
  }

  for (name in names(JASPAR_PWM)){
    # Combine motifs of JASPAR20202 and chromVARmotif
    motifs[name] <- JASPAR_PWM[name]
  }
  tf2motifs <- data.frame(motif = character(),
  # Initiate final TF motifs table
                          tf = character(),
                          stringsAsFactors = FALSE)
  for (i in seq_along(names(motifs))){                   # Fill TF motif table
    tfs <- strsplit(names(motifs)[i], "::")[[1]]
    # splitting TFs that are given as "name1::name2"
    for (tf in tfs){
      tf <- strsplit(tf, "(", fixed = TRUE)[[1]][1]
      # only keeping <NAME> in identifier "<NAME>(var.n)"
      tf2motifs <- rbind(tf2motifs, data.frame(motif = names(motifs)[i], tf = tf))
    }
  }

  return(new("motifs_db",
             tf2motifs = tf2motifs,
             motifs = motifs))
}
