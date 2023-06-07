## ----setup, include=FALSE-----------------------------------------------------
#knitr::opts_chunk$set(eval = FALSE)
devtools::install_github("cantinilab/HuMMuS")


## -----------------------------------------------------------------------------
library(HuMMuS)
library(reticulate)


## -----------------------------------------------------------------------------
# Load the Chen dataset, which is a Seurat object containing scRNA-seq and scATAC-seq data
data("chen_dataset_subset")
chen_dataset_subset

# Create an hummus object from seurat object
hummus <- as(chen_dataset_subset, 'hummus_object')

a <- chen_dataset


## ----warning=FALSE------------------------------------------------------------
# get human genome annotation from EndDb data
genome_annotations <- get_genome_annotations(
  ensdb_annotations = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)


## -----------------------------------------------------------------------------
# We add annotation to the ATAC/peak assay
# (will be used to define peak-gene links)
Signac::Annotation(hummus@assays$peaks) <- genome_annotations

# Load TF motifs from JASPAR2020 and chromVARmotifs in hummus object
hummus@motifs_db <- get_tf2motifs() # by default human motifs


## -----------------------------------------------------------------------------
hummus <- bipartite_tfs2peaks(
              hummus_object = hummus,
              tf_expr_assay = "RNA", # use to filter TF on only expressed TFs,
                                     # if NULL, all TFs with motifs are used
              peak_assay = "peaks",
              tf_multiplex_name = "TF",
              genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
              )


## -----------------------------------------------------------------------------
hummus <- bipartite_peaks2genes(
                      hummus_object = hummus,
                      gene_assay = "RNA",
                      peak_assay = "peaks",
                      store_network = FALSE,
                      )


## -----------------------------------------------------------------------------
hummus <- compute_tf_network(hummus,
                            gene_assay = "RNA", # default = None ;
                                                #If a assay is provided,
                                                # only the TFs that are present
                                                # will be considered
                            verbose = 1,
                            multiplex_name = "TF",
                            tf_network_name = "TF_network")


## -----------------------------------------------------------------------------
hummus <- compute_gene_network(
              hummus,
              gene_assay = "RNA",
              method = "GENIE3",
              verbose = 1,
              number_cores = 5, # GENIE3 method can be ran
                                # parallelised on multiple cores
              store_network = FALSE, # by default : FALSE
                                     # each network can be saved 
                                     # when computed with hummus
              output_file = "gene_network.tsv")


## -----------------------------------------------------------------------------
hummus <- compute_atac_peak_network(hummus,
              atac_assay = "peaks",
              verbose = 1,
              genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
              store_network = FALSE)  


## -----------------------------------------------------------------------------
save_multilayer(hummus = hummus, "chen_multilayer")


## -----------------------------------------------------------------------------
grn <- define_grn(
  hummus,
  multilayer_f = "chen_multilayer",
  njobs = 5
  )

## -----------------------------------------------------------------------------
head(grn)


## -----------------------------------------------------------------------------
enhancers <- define_enhancers(
  hummus,
  multilayer_f = "chen_multilayer",
  njobs = 5
  )

## -----------------------------------------------------------------------------
head(enhancers)


## -----------------------------------------------------------------------------
binding_regions <- define_binding_regions(
  hummus,
  multilayer_f = "chen_multilayer",
  njobs = 5
  )

## -----------------------------------------------------------------------------
head(binding_regions)


## -----------------------------------------------------------------------------
target_genes <- define_target_genes(
  hummus,
  multilayer_f = "chen_multilayer",
  njobs = 5
  )

## -----------------------------------------------------------------------------
head(target_genes)

