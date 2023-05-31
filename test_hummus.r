devtools::load_all("../hummus_package")
library(reticulate)
use_condaenv("base")


multilayer_folder <- "examples/multilayer_1"
   #### Create hummus object ####
# Data
scRNA <- read.table("data/real_example/hESC_Chen_scRNA.tsv")
scRNA <- scRNA[1:700, ]
scATAC <- read.table("data/real_example/hESC_Chen_scATAC_bin.tsv")
scATAC <- scATAC[1:800, ]
# Create seurat object with scRNA-seq data
seurat <- SeuratObject::CreateSeuratObject(scRNA)
# Add scATAC-seq data to seurat object
seurat[["peaks"]] <- Signac::CreateChromatinAssay(scATAC, sep = c(":", "-"))


###############################################
# Create hummus object from seurat object     #
###############################################
aa <- Sys.time()

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
hummus <- bipartite_peaks2genes(
                      hummus_object = hummus,
                      gene_assay = "RNA",
                      peak_assay = "peaks",
                      store_network = FALSE,
                      )

# Add bipartite between tfs and peaks data to hummus object
hummus <- bipartite_tfs2peaks(
              hummus_object = hummus,
              tf_expr_assay = "RNA", #use to filter TF on only expressed TFs, if NULL, all TFs with motifs are used
              peak_assay = "peaks",
              tf_network_name = "TF",
              genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
              )

# Compute layer networks and add it to hummus object
hummus <- compute_tf_network(hummus,
                            #gene_assay = "RNA",
                            verbose = 1,
                            multiplex_name = "TF",
                            tf_network_name = "TF_network",
                            )


# Compute and save gene network from scRNA-seq w/ GENIE3
hummus <- compute_gene_network(hummus,
                              gene_assay = "RNA",
                              method = "GENIE3",
                              verbose = 1,
                              number_cores = 8,
                              store_network = FALSE)

hummus <- compute_atac_peak_network(hummus,
                                    atac_assay = "peaks",
                                    verbose = 1,
        genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
                                    store_network = FALSE)

# Save the mulilayer in a classic hierarchical structure
save_multilayer(hummus = hummus, multilayer_folder)

cat("Total multilayer time construction : ", Sys.time() - aa, "\n")
bb <- Sys.time() - aa

# Define GRN
grn <- define_grn(
  hummus,
  multilayer_f = multilayer_folder,
  njobs = 5
  )


#"A" :
  #  "bipartite":
    #  "atac_rna.tsv"
    #  "tf_peaks.tsv"
  #  "multiplex":
    #  "gene_network.tsv"
    #  "tf_network.tsv"
    #  "atac_network.tsv"
  #  "seeds"
  #  "config"