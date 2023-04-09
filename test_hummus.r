devtools::load_all("../hummus_package")

# Data
scRNA <- read.table("data/real_example/hESC_Chen_scRNA.tsv")
scRNA <- scRNA[1:700, ]

scATAC <- read.table("data/real_example/hESC_Chen_scATAC_bin.tsv")
scATAC <- scATAC[1:1000, ]

# Create seurat object with scRNA-seq data
seurat <- SeuratObject::CreateSeuratObject(scRNA)
# Add scATAC-seq data to seurat object
seurat[["peaks"]] <- Signac::CreateChromatinAssay(scATAC, sep = c(":", "-"))

# Fetch genome annotations
#genome_annotations <- get_genome_annotations(
 # ensdb_annotations = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)

# Add genome annotations to seurat object
Signac::Annotation(seurat@assays$peaks) <- genome_annotations

# Create hummus object
hummus <- hummus_object(seurat)

# Load TF motifs from JASPAR2020 and chromVARmotifs
hummus@motifs_db <- get_tf2motifs()

# Compute and save peak-gene links
hummus@multilayer@bipartites$atac_rna <- bipartite_peaks2genes(
                      seurat_object = seurat,
                      gene_assay = "RNA",
                      peak_assay = "peaks",
                      store_bipartite = FALSE,
                      )



# Get TFs expressed in scRNA-seq data and having known binding motifs
tfs <- get_tfs(hummus, assay = "RNA", store_tfs = FALSE)

hummus@multilayer@bipartites$tf_peaks <- bipartite_tfs2peaks(
                      hummus_object = hummus,
                      tf_gene_assay = "RNA",
                      peak_assay = "peaks",
                      genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
                      store_bipartite = FALSE,
                      )


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