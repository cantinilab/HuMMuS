#' Title
#'
#' @param l_genes
#' @param l_peaks
#' @param peak_sep1
#' @param peak_sep2
#' @param store_bipartite
#' @param output_file
#' @param genome
#' @param gene.range
#' @param motifs
#' @param tf2motifs
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
Bipartite_TFs2Peaks <- function(
  l_genes,
  l_peaks,
  peak_sep1=":",
  peak_sep2="-",
  store_bipartite=TRUE,
  output_file,
  genome,
  gene.range,
  motifs,
  tf2motifs,
  verbose=1){

  # Build up object to determine TF-peak links and peak-gene links
  rna  = data.frame(features=l_genes)                                           # List of genes present in scRNA
  rna[,c("dummy_cell1", "dummy_cell2")] <- 1
  rownames(rna) <- rna$features
  rna$features <- NULL
  atac  = data.frame(features=l_peaks)                                          # List of peaks present in scATAC
  atac[,c("dummy_cell1", "dummy_cell2")] <- 1
  rownames(atac) <- atac$features
  atac$features <- NULL
  seurat = CreateSeuratObject(rna)                                              # Create seurat object with genes
  seurat[['peaks']] <- CreateChromatinAssay(atac, sep=c(peak_sep1,peak_sep2))   # Combine genes and peaks in a seurat object
  Annotation(seurat@assays$peaks) <- gene.range# Add genome annotations to seurat object

  #cand_ranges <- object@grn@regions@ranges
  motif_pos <- Signac::AddMotifs(
    object = seurat[['peaks']],
    genome = genome,
    pfm = motifs,  #add verbose options
  )
  #regons_motifs <- AddMotifs(seurat[['peaks']], genome, motifs)

  #### The 17 following lines are inspired from the Pando package : https://github.com/quadbiolab/Pando/blob/main/R/regions.R
  # Add TF info for motifs
  if (verbose>0){
    print('Adding TF info', verbose=verbose)
  }
  if (!is.null(tf2motifs)){
    tf2motifs <- tibble(tf2motifs)
  } else {
    utils::data(tf2motifs, envir = environment())
  }

  # Spread dataframe to sparse matrix
  tf2motifs <- tf2motifs %>% dplyr::select('motif'=1,'tf'=2) %>%
    distinct() %>% mutate(val=1) %>%
    tidyr::pivot_wider(names_from = 'tf', values_from=val, values_fill=0) %>%
    column_to_rownames('motif') %>% as.matrix() %>% Matrix::Matrix(sparse=TRUE)
  tfs_use <- intersect(rownames(GetAssay(seurat, 'RNA')), colnames(tf2motifs))
  if (length(tfs_use)==0){
    stop('None of the provided TFs were found in the dataset. Consider providing a custom motif-to-TF map as `motif_tfs`')
  }

  TFs_Peaks = motif_pos@motifs@data %*% tf2motifs[, tfs_use] # Get TF peak links
  TFs_Peaks = TFs_Peaks[,colnames(TFs_Peaks) %in% l_genes]                      # Keep only the TFs that are in our scRNA-seq dataset
  l_tfs2peaks <- expand.grid(rownames(TFs_Peaks),colnames(TFs_Peaks))[as.vector(TFs_Peaks>0),] # TF-peak links
  colnames(l_tfs2peaks) <- c("peak","TF")                                                      # set column names

  if(store_bipartite){
    write.table(l_tfs2peaks, output_file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  }

  return(l_tfs2peaks)
}




#' Title
#'
#' @param genes
#' @param peaks
#' @param peak_sep1
#' @param peak_sep2
#' @param gene.range
#' @param peak_to_gene_method
#' @param upstream
#' @param downstream
#' @param only_tss
#' @param store_bipartite
#' @param output_file
#' @param peaks2genes_filename
#'
#' @return
#' @export
#'
#' @examples
Bipartite_Peaks2Genes <- function(
  genes,
  peaks,
  peak_sep1,
  peak_sep2,
  gene.range,
  peak_to_gene_method='Signac',
  upstream=500,
  downstream=500,
  only_tss=TRUE,
  store_bipartite=TRUE,
  output_file,
  peaks2genes_filename){

  # Build up object to determine peak-gene links
  rna  = data.frame(features=genes)                                           # List of genes present in scRNA
  rna[,c("dummy_cell1", "dummy_cell2")] <- 1
  rownames(rna) <- rna$features
  rna$features <- NULL
  atac  = data.frame(features=peaks)                                          # List of peaks present in scATAC
  atac[,c("dummy_cell1", "dummy_cell2")] <- 1
  rownames(atac) <- atac$features
  atac$features <- NULL
  seurat = CreateSeuratObject(rna)                                              # Create seurat object with genes
  seurat[['peaks']] <- CreateChromatinAssay(atac, sep=c(peak_sep1,peak_sep2))   # Combine genes and peaks in a seurat object
  Annotation(seurat@assays$peaks) <- gene.range# Add genome annotations to seurat object


    #Params <- Params(seuratPlus)                                                  # Get variables (name of peak_assay and rna_assay and boolean exlude_exon)
  regions <- NetworkRegions(seuratPlus)                                         # Get network regions
  gene_annot <- Signac::Annotation(seuratPlus[['peaks']])             # gene annotation

  peaks_near_gene <- find_peaks_near_genes(peaks = seurat[['peaks']]@ranges,              # Find candidate regions near gene bodies
                                           method = peak_to_gene_method,
                                           genes = gene_annot,
                                           upstream = upstream,
                                           downstream = downstream,
                                           only_tss = only_tss)

  peaks2gene <- aggregate_matrix(t(peaks_near_gene),                            # Get gene to peak matrix
                                 groups=colnames(peaks_near_gene),
                                 fun="sum")
  peaks2gene = peaks2gene[rownames(peaks2gene) %in% genes,]                   # Keep only the genes that are in our scRNA-seq dataset
  peaks2gene = peaks2gene[rowSums(peaks2gene)!=0, colSums(peaks2gene)!=0]       # Remove rows/cols with only zeros
  l_peaks2genes <- expand.grid(rownames(peaks2gene),colnames(peaks2gene))[as.vector(peaks2gene>0),] # peak-gene links
  colnames(l_peaks2genes) <- c("gene","peak")                                                      # set column names

  if(store_bipartite){
    write.table(l_peaks2genes, output_file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  }

  # Return list the two edgelists containing TF-peak and peak-gene links
  return(l_peaks2genes)
}
