devtools::load_all("../hummus_package")
   #### Create hummus object ####
# Data
scRNA <- read.table("data/real_example/hESC_Chen_scRNA.tsv")
scRNA <- scRNA[1:700, ]
scATAC <- read.table("data/real_example/hESC_Chen_scATAC_bin.tsv")
scATAC <- scATAC[1:1000, ]
# Create seurat object with scRNA-seq data
seurat <- SeuratObject::CreateSeuratObject(scRNA)
# Add scATAC-seq data to seurat object
seurat[["peaks"]] <- Signac::CreateChromatinAssay(scATAC, sep = c(":", "-"))


###############################################
# Create hummus object from seurat object     #
###############################################
hummus <- hummus_object(seurat)

   #### Add annotations to hummus object
# Fetch genome annotations
genome_annotations <- get_genome_annotations(
  ensdb_annotations = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
# Add genome annotations to seurat object
Signac::Annotation(hummus@assays$peaks) <- genome_annotations
# Load TF motifs from JASPAR2020 and chromVARmotifs in hummus object
hummus@motifs_db <- get_tf2motifs()

   #### Compute bipartites and add it to hummus object
# Add bipartite between rna-genes and peaks data to hummus object
hummus@multilayer@bipartites$atac_rna <- bipartite_peaks2genes(
                      seurat_object = hummus,
                      gene_assay = "RNA",
                      peak_assay = "peaks",
                      store_bipartite = FALSE,
                      )

# Add bipartite between tfs and peaks data to hummus object
hummus@multilayer@bipartites$tf_peaks <- bipartite_tfs2peaks(
              hummus_object = hummus,
              tf_expr_assay = "RNA", #use to filter TF on only expressed TFs, if NULL, all TFs with motifs are used
              peak_assay = "peaks",
              tf_network_name = "TF_network",
              genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
              )

   #### Compute layer networks and add it to hummus object
hummus <- compute_tf_network(hummus,
                            gene_assay = "RNA",
                            verbose = 1,
                            multiplex_name = "TF",
                            )


# Compute and save gene network from scRNA-seq w/ GENIE3
hummus <- compute_gene_network(hummus,
                              gene_assay = "RNA",
                              method = "GENIE3",
                              verbose = 1,
                              number_cores = 5,
                              store_network = TRUE,
                              output_file = "gene_network.tsv")


hummus <- compute_gene_network(hummus,
                              gene_assay = "RNA",
                              method = "GENIE3",
                              verbose = 1,
                              number_cores = 5,
                              multiplex_name = "peak_network",
                              )

# Save the mulilayer in a classic hierarchical structure
save_multilayer(hummus = hummus, "a")
#"A" :
  #  "bipartites":
    #    "atac_rna"
    #    "tf_peaks"
  #  "multiplex":
    #    "gene_network"
    #    #"tf_network"
    #    #"atac_network"
  #  "seeds"
  #  "config"









# Compute and save TF network
#compute_TF_network(
 #     tfs = tfs,
  #    output_file = "data/real_example/multilayer/multiplex/hESC_Chen_TF_network.tsv")

# Compute and save gene network
#compute_gene_network(
 #     scRNA = scRNA,
  #    tfs = tfs,
   #   method = "GENIE3",
    #  output_file = "data/real_example/multilayer/multiplex/hESC_Chen_gene_network.tsv",
     # verbose = 1,
      #number_cores = 5)

# Compute and save peak network
#compute_atac_peak_network(
 #     scATAC = scATAC,
  #    output_file = "data/real_example/multilayer/multiplex/hESC_Chen_ATAC_network.tsv",
   #   genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    #  verbose = 1)