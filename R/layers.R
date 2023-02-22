############################################################################################################################################################
## compute TF network
############################################################################################################################################################

#' Title
#'
#' @param organism
#' @param list_TF
#' @param store_network
#' @param output_file
#' @param source_target
#'
#' @return
#' @export
#'
#' @examples
compute_TF_network <- function(organism=9606,
                               list_TF=NA, store_network=TRUE, output_file, source_target='AND'){

   TF_PPI <- import_post_translational_interactions(
    organism = organism, partners = list_TF, source_target = source_target
  )

  Network_TF = TF_PPI[,c(3,4)]

  if (store_network==TRUE){
    write.table(Network_TF, output_file,                                                  # Store edgelist
                col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  }

  return(Network_TF)
}

############################################################################################################################################################
## compute gene network
############################################################################################################################################################

#' Title
#'
#' @param scRNA
#' @param TFs
#' @param method
#' @param store_network
#' @param output_file
#' @param threshold
#' @param number_cores
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
compute_gene_network <- function(scRNA, TFs, method="GENIE3", store_network=TRUE, output_file, threshold=0.0, number_cores=1, verbose=1){
  if(method=="GENIE3"){
    a = Sys.time()
    weightMat = GENIE3(as.matrix(scRNA), regulators = TFs, nCores=number_cores)   # Infer network

    if (verbose>0){
      print(paste("Gene network construction time:",Sys.time()-a))
    }

    linkList = getLinkList(weightMat)     # Get edge list
    gene_network = linkList[which(linkList$weight>threshold),]

    if (store_network==TRUE){
      write.table(gene_network, output_file,                                                           # Store edgelist
                  col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    }

    return(gene_network)                                     # Return edges with weight > 0
  }
  else{
    return("no gene network computed!")
  }
}

#' Compute netwok from atac data
#'
#' This function will create a network from atac data (or in theory any data
#' wtih peaks coordinates as features).
#' Different method should be implemented at some point (e.g. RENIN),
#' for now Cicero is still the reference and only method available
#'
#' Method descriptions :
#'  1. Cicero
#'      Use patial corelation between peaks that are in a given window (e.g. :
#'      less distant than 500K base pairs)
#'
#' @param scATAC
#' @param genome
#' @param method
#' @param store_network
#' @param output_file
#' @param threshold
#' @param number_cluster
#' @param number_sample
#' @param seed
#' @param verbose
#' @param window
#'
#' @return
#' @export
#'
#' @examples
compute_atac_peak_network <- function(
    scATAC,
    genome=BSgenome.Hsapiens.UCSC.hg38,
    method="cicero",
    store_network=TRUE,
    output_file,
    threshold=0.0,
    number_cluster=50,
    number_sample=100,
    seed=2021,
    verbose=1,
    window = 5e+05
    ){

  if(method=="cicero"){
    chromosome_sizes <- data.frame(V1=genome@seqinfo@seqnames,    # obtain chromosome sizes
                                   V2=genome@seqinfo@seqlengths)
    # Matrix to edgelist
    ACC <- reshape2::melt(as.matrix(scATAC))
    colnames(ACC) <- c("V1","V2","V3")
    ACC$V1 <- gsub("_", "-", ACC$V1)
    # Run cicero
    input_cds <- make_atac_cds(ACC, binarize = TRUE)                         # Create CDS object
    set.seed(seed)
    print('test')
    if(length(which(colSums(exprs(input_cds))==0))==0)                       # It is required that there is no empty cell
    {
      input_cds <- estimate_size_factors(input_cds)                            # Calculating size factors using default method = mean-geometric-mean-total
      input_cds <- preprocess_cds(input_cds, method = "LSI")                   # Preprocessing using LSI
      input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP',      # Dimensionality reduction
                                    preprocess_method = "LSI")
    }
    else{
      print("Error: there is at least one cell with no signal.")
    }
    print(input_cds)
    umap_coords <- reducedDims(input_cds)$UMAP                                # Get reduced (UMAP) coordinates
    cicero_cds <- make_cicero_cds(input_cds,                                  # Create a Cicero CDS object
                                  reduced_coordinates = umap_coords,
                                  k = number_cluster,  #number neighbors         # Default = 50
                                  summary_stats = NULL,         # Default
                                  size_factor_normalize = TRUE, # Default
                                  silent = FALSE)               # Default
    print('test 2')
    a = Sys.time()
    cicero <- run_cicero(cds = cicero_cds,                                   # Infer peak-links
                         genomic_coords = chromosome_sizes, # mouse.mm10.genome,
                         window = window,             # Default
                         silent = FALSE,             # Default
                         sample_num = number_sample) # Default = 100
    print('test 3')
    if (verbose<0){
      print(paste("Peak network construction time:",Sys.time()-a))}
    if(store_network){
      write.table(cicero, cicero_nw_filename, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    }
    # Remove NAs, double edges, and edges with weight <=0
    if(length(which(is.na(cicero$coaccess)))>threshold)                                # Check for coaccess = NA
    {
      cicero <- cicero[which(!is.na(cicero$coaccess)),]                        # Remove NAs
    }
    cicero$temp <- NA                                                          # Helper column to check and remove double edges
    my_cols <- which(as.character(cicero$Peak1) <= as.character(cicero$Peak2))
    cicero$temp[my_cols] <- paste(cicero$Peak1[my_cols],cicero$Peak2[my_cols],sep = ";")
    my_cols <- which(as.character(cicero$Peak1) > as.character(cicero$Peak2))
    cicero$temp[my_cols] <- paste(cicero$Peak2[my_cols],cicero$Peak1[my_cols],sep = ";")
    cicero <- cicero[with(cicero, order(temp, decreasing = TRUE)),]            # Sort table according to temp-column (each entry appears double)
    rownames(cicero) <- c(1:dim(cicero)[1])
    A <- as.character(cicero$Peak1[seq(1,dim(cicero)[1],2)])
    Anum <- round(cicero$coaccess[seq(1,dim(cicero)[1],2)],10)
    B <- as.character(cicero$Peak2[seq(2,dim(cicero)[1],2)])
    Bnum <- round(cicero$coaccess[seq(2,dim(cicero)[1],2)],10)
    #length(which(A==B & Anum==Bnum))                                          # Each edge appears twice with same coaccess score (rounded to 10 digits after comma)
    cicero <- cicero[seq(1,dim(cicero)[1],2),]                                 # Remove double edges
    cicero$temp <- NULL                                                        # Remove helper column
    cicero <- cicero[with(cicero, order(coaccess, decreasing = TRUE)),]        # Sort
    rownames(cicero) <- c(1:dim(cicero)[1])
    cicero$Peak1 <- gsub("_","-",cicero$Peak1)                                 # Peak names 2x"-" to match bipartites
    cicero$Peak2 <- gsub("_","-",cicero$Peak2)                                 # Peak names 2x"-" to match bipartites########? 2x"-" or 2x"_"

    if (verbose<0){
      print(paste(dim(Peak_Network)[1], "peak edges with weight > ", as.string(threshold)))}

    peak_network = cicero[which(cicero$coaccess>threshold),]

    if(store_network){
      write.table(peak_network, outputfilename,
                  col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t") # Store edgelist
    }

    # Return peak network including edges with positive coaccess score
    print('test')
    return(peak_network)
  }
  else{
    return("no peak network computed!")
  }
}
