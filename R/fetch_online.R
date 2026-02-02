#'  Fetch online genome annotations from Ensembldb database
#'
#' @param EnsDb_annotations (EndsDb object) - Ensembldb database (default: EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
#'
#' @return gene_range (GRanges object) - Genome annotations
#' @export
#'
#' @examples gene_range = get_genome_annotations(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
get_genome_annotations <- function(
  ensdb_annotations = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
  ) {
  # Get genome annotations from Ensembldb database
  gene_range <- Signac::GetGRangesFromEnsDb(ensdb_annotations)

  ucsc.levels <- stringr::str_replace(
    string = paste("chr", Signac::seqlevels(gene_range), sep = ""),
    pattern = "chrMT",
    replacement = "chrM") # Change chromosome names to UCSC format

  Signac::seqlevels(gene_range) <- ucsc.levels
  # check if Signac is the good package

  return(gene_range) # Return genome annotations
}

#' Fetch online TF motifs from JASPAR2020 and chromVARmotifs
#'
#' @param species (character) - Species name (default: "human")
#'
#' @return motifs_db (motifs_db object) - TF2motifs + motifs PWMs
#' @export
#'
#' @examples motifs_db = get_tf2motifs(species = "human")
get_tf2motifs <- function(species = "human", download_folder = NULL ) {
  if( ! is.null(download_folder == NULL) ){
    download_folder = getwd()
  }
  
  #TF motifs using the union of databases: JASPAR and cis-BP
  # included in chromVAR
  getMatrixSet <- TFBSTools::getMatrixSet

  # If species is human or mouse
  if (species == "human") {
    # Parameters for JASPAR2020
    opts <- list(collection = "CORE",
              species    = "Homo sapiens",
              all_versions = FALSE)
    JASPAR_PWM <- TFBSTools::toPWM(getMatrixSet(JASPAR2020::JASPAR2020, opts))
    # Load data from JASPAR2020
    # Load data from chromVARmotifs
    
    # Original data accessible at https://github.com/GreenleafLab/chromVARmotifs
    url = "https://github.com/GreenleafLab/chromVARmotifs/raw/refs/heads/master/data/human_pwms_v2.rda"
    dpath = file.path( download_folder, 'human_pwms_v2.rda' )
    download.file(url, dpath)
    load( file=dpath )
    
    # Load data from chromVARmotifs
    motifs <- human_pwms_v2
    # Motifs from chromVARmotifs
  } else if (species == "mouse") {
    # Parameters for JASPAR2020
    opts <- list(collection = "CORE",
              species    = "Mus musculus",
              all_versions = FALSE)
    JASPAR_PWM <- TFBSTools::toPWM(getMatrixSet(JASPAR2020::JASPAR2020, opts))
    # Load data from JASPAR2020
    url = "https://github.com/GreenleafLab/chromVARmotifs/raw/refs/heads/master/data/mouse_pwms_v2.rda"
    dpath = file.path( download_folder, 'mouse_pwms_v2.rda' )
    download.file(url, dpath)
    load( file=dpath )
    
# Load data from chromVARmotifs
    # Original data accessible at https://github.com/GreenleafLab/chromVARmotifs
    motifs <- mouse_pwms_v2
    # Motifs from chromVARmotifs
  }

  for (name in names(JASPAR_PWM)){
    # Combine motifs of JASPAR20202 and chromVARmotif
    motifs[name] <- JASPAR_PWM[name]
  }

  # Initiate final TF motifs table
  tf2motifs <- data.frame(motif = character(),
                          tf = character(),
                          stringsAsFactors = FALSE)
  for (i in seq_along(TFBSTools::name(motifs))){  # Fill TF motif table
  # TFBSTools::name(motifs) returns names of TFs associated to each PWMatrix
    tfs <- strsplit(TFBSTools::name(motifs)[i], "::")[[1]]
    # splitting TFs that are given as "name1::name2"
    for (tf in tfs){
      tf <- strsplit(tf, "(", fixed = TRUE)[[1]][1]
      # only keeping <NAME> in identifier "<NAME>(var.n)"
      tf2motifs <- rbind(tf2motifs, data.frame(motif = names(motifs)[i],
                                                             tf = tf))
    }
  }

  return(new("motifs_db",
             tf2motifs = tf2motifs,
             motifs = motifs,
             tfs = unique(tf2motifs$tf))) # Return motifs_db <- TF2motifs + motifs PWMs
}
