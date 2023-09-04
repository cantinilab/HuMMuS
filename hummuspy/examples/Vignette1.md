Predict GRN from scRNA+scATAC data (Chen 2018 dataset)
================
Trimbour Rémi
2023-05-16

## General description of the vignette

##### 0. Preparation of the data

##### 1. Initialisation of a hummus object

##### 2. Construction of the multilayer

We want to link TFs to peaks (atac regions), then peaks (atac regions)
to genes. We thus need 2 bipartites:

- 2.1. **TF - peaks bipartite**
- 2.2. **peaks - genes bipartite**

We also need one layer per type of features : TF, peaks and genes.

- 2.3. **TFs/Proteins layer**
- 2.4 **Peaks/ATAC layer**
- 2.5 **Genes/RNA layer**

##### 3. Analyse multilayer and define gene regulatory network (GRN)\$

This part requires first have hummuspy python package (pip install
hummuspy) installed, and second to save locally the multilayer files

- 3.1. Save the multilayer
- 3.2. Create GRN output

## Download the single-cell data

The data used in this tutorial can be downloaded at : ibens repo link

## 0. Import single-cell data

``` r
scRNA <- read.table("data/real_example/hESC_Chen_scRNA.tsv")
scATAC <- read.table("data/real_example/hESC_Chen_scATAC_bin.tsv")
scRNA <- scRNA[1:700, ]
scATAC <- scATAC[1:1000, ]
```

## 1. Initialisation of HuMMuS object

HuMMuS R objects are instances developed on top of seurat objects. It
means it’s created from a seurat object and the contained assays can be
accessed the same way.

Additionally, it contains a motifs_db object, providing tf motifs
informations, and a multilayer objects, that will be completed while
going through this tutorial. It will mostly include :

- list of multiplex networks (one per modality)
- list of bipartites (one per connection between layers)

#### 1.1. Transform data into a hummus object

``` r
# Create seurat object with scRNA-seq data
seurat <- SeuratObject::CreateSeuratObject(scRNA)

# Add scATAC-seq data to seurat object
seurat[["peaks"]] <- Signac::CreateChromatinAssay(scATAC, sep = c(":", "-"))

# Create an hummus object from seurat object
hummus <- hummus_object(seurat)
```

#### 1.2. Add genome annotations to hummus object

Fetch genome annotations

``` r
# get human genome annotation from EndDb data
genome_annotations <- get_genome_annotations(
  ensdb_annotations = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
```

Add genome annotations to hummus/seurat object

``` r
# We add annotation to the ATAC/peak assay
# (will be used to define peak-gene links)
Signac::Annotation(hummus@assays$peaks) <- genome_annotations
# Load TF motifs from JASPAR2020 and chromVARmotifs in hummus object
hummus@motifs_db <- get_tf2motifs() # by default human motifs
```

## 2. Construction of the multilayer

### Compute bipartites and add it to hummus object

#### 2.1. TF - peaks bipartite reconstruction

``` r
hummus <- bipartite_tfs2peaks(
              hummus_object = hummus,
              tf_expr_assay = "RNA", # use to filter TF on only expressed TFs,
                                     # if NULL, all TFs with motifs are used
              peak_assay = "peaks",
              tf_multiplex_name = "TF",
              genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
              )
```

    ## Computing TF-peak bipartite
    ##   13 TFs expressed

    ## Building motif matrix

    ## Finding motif positions

    ## Creating Motif object

    ##  Adding TF info
    ##  Returning TF-peak links as bipartite object
    ## no peak layer name provided, using peak_assay name

#### 2.2. Genes - peaks bipartite reconstruction

``` r
hummus <- bipartite_peaks2genes(
                      hummus_object = hummus,
                      gene_assay = "RNA",
                      peak_assay = "peaks",
                      store_network = FALSE,
                      )
```

### Compute layer networks and add it to hummus object

To illustrate the possible level of parameters specifications, each
network will be reconstruct with different personnalisation levels.

#### 2.3. Compute the TF network from OmniPath database

``` r
hummus <- compute_tf_network(hummus,
                            gene_assay = "RNA", # default = None ;
                                                #If a assay is provided,
                                                # only the TFs that are present
                                                # will be considered
                            verbose = 1,
                            multiplex_name = "TF",
                            tf_network_name = "TF_network")
```

    ## Computing TF network...
    ##   13 TFs expressed
    ##  TF network construction time: 0.9180498 
    ## No TF-TF edges from Omnipath for the given parameters.
    ## 
    ##         You can try to change the source_target parameter to 'OR' to get
    ##         TF-other protein interactions.
    ##  Or try to import a network 
    ##         computed externally.
    ##  Right now, a network with all TFs connected
    ##         to a fake node is created, for HuMMuS analysis.
    ##  It has no biological
    ##         meaning but will allow to run the pipeline as if no edges were present.
    ##         
    ##  Creating new multiplex :  TF

#### 2.4. Compute gene network from scRNA-seq w/ GENIE3

``` r
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
```

    ## Computing gene network with  GENIE3  ...
    ##  No TFs list provided, fetching from hummus object...
    ##   13 TFs expressed
    ##  Gene network construction time: 1.134216 
    ##  Creating new multiplex :  RNA

#### 2.5. Compute the peak network from scATAC-seq w/ Cicero

``` r
hummus <- compute_atac_peak_network(hummus,
              atac_assay = "peaks",
              verbose = 1,
              genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
              store_network = FALSE)  
```

    ## Overlap QC metrics:
    ## Cells per bin: 50
    ## Maximum shared cells bin-bin: 44
    ## Mean shared cells bin-bin: 7.04584942084942
    ## Median shared cells bin-bin: 0

    ## Warning in cicero::make_cicero_cds(input_cds, reduced_coordinates =
    ## umap_coords, : On average, more than 10% of cells are shared between paired
    ## bins.

    ## [1] "Starting Cicero"
    ## [1] "Calculating distance_parameter value"

    ## Warning in estimate_distance_parameter(cds, window = window, maxit = 100, :
    ## Could not calculate sample_num distance_parameters (2 were calculated) - see
    ## documentation details

    ## [1] "Running models"
    ## [1] "Assembling connections"
    ## [1] "Successful cicero models:  226"
    ## [1] "Other models: "
    ## 
    ## Zero or one element in range 
    ##                        13215 
    ## [1] "Models with errors:  0"
    ## [1] "Done"
    ## 
    ##  71 peak edges with a coaccess score > 0 were found.
    ## Peak network construction time: 36.75066 Creating new multiplex :  peaks

## Analyse of the multilayer and definition of GRN

#### 3.1. Save the mulilayer in a classic hierarchical structure

``` r
save_multilayer(hummus = hummus, "chen_multilayer")
```

    ## Warning in dir.create(folder_name): 'chen_multilayer' already exists

    ## Warning in dir.create(paste0(folder_name, "/", multiplex_folder)):
    ## 'chen_multilayer/multiplex' already exists

    ## Warning in dir.create(paste0(folder_name, "/", bipartite_folder)):
    ## 'chen_multilayer/bipartite' already exists

    ## Warning in dir.create(paste0(folder_name, "/", seed_folder)):
    ## 'chen_multilayer/seed' already exists

    ## Warning in dir.create(paste0(folder_name, "/", config_folder)):
    ## 'chen_multilayer/config' already exists

    ## Warning in dir.create(paste0(folder_name, "/", multiplex_folder, "/",
    ## multiplex_name)): 'chen_multilayer/multiplex/TF' already exists

    ## Multiplex of  1  networks with 14 features.
    ##  Networks names:  TF_network[1] "TF TF_network"

    ## Warning in dir.create(paste0(folder_name, "/", multiplex_folder, "/",
    ## multiplex_name)): 'chen_multilayer/multiplex/RNA' already exists

    ## Multiplex of  1  networks with 700 features.
    ##  Networks names:  RNA_GENIE3[1] "RNA RNA_GENIE3"

    ## Warning in dir.create(paste0(folder_name, "/", multiplex_folder, "/",
    ## multiplex_name)): 'chen_multilayer/multiplex/peaks' already exists

    ## Multiplex of  1  networks with 116 features.
    ##  Networks names:  peak_network_cicero[1] "peaks peak_network_cicero"

#### 3.2. Define GRN

``` r
grn <- define_grn(
  hummus,
  multilayer_f = "chen_multilayer",
  njobs = 5
  )
```

``` r
head(grn)
```

    ##   layer                  path_layer        score  gene        tf
    ## 1    TF multiplex/TF/TF_network.tsv 0.0219787526 ABHD8    ARID3B
    ## 2    TF multiplex/TF/TF_network.tsv 0.0047128575 ABHD8 fake_node
    ## 3    TF multiplex/TF/TF_network.tsv 0.0002252904 ABHD8    ARID3A
    ## 4    TF multiplex/TF/TF_network.tsv 0.0001813984 ABHD8     BACH1
    ## 5    TF multiplex/TF/TF_network.tsv 0.0001541393 ABHD8       BBX
    ## 6    TF multiplex/TF/TF_network.tsv 0.0001238019 ABHD8      ATF4
