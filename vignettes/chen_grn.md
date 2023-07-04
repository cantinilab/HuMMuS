Predict GRN from scRNA+scATAC data (Chen 2018 dataset)
================
Trimbour Rémi
2023-05-16

## General description of the vignette

##### 0. Preparation of the environment

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

##### 3. Analyse multilayer and define gene regulatory network (GRN)

This part requires first have hummuspy python package (pip install
hummuspy) installed, and second to save locally the multilayer files

- 3.1. Save the network files
- 3.2. Create GRN output
- 3.3. Create enhancers output
- 3.4. Create binding regions output
- 3.5. Create target genes output

## Download the single-cell data

The data used in this tutorial can be [downloaded
here](https://figshare.com/account/home#/projects/168899)

## 0. Setting up the environment

``` r
library(HuMMuS)
```

    ## The legacy packages maptools, rgdal, and rgeos, underpinning this package
    ## will retire shortly. Please refer to R-spatial evolution reports on
    ## https://r-spatial.org/r/2023/05/15/evolution4.html for details.
    ## This package is now running under evolution status 0

    ## 

``` r
library(reticulate)
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

### 1.1. Transform data into a hummus object

``` r
# Load the Chen dataset, which is a Seurat object containing scRNA-seq and scATAC-seq data
data("chen_dataset_subset")
chen_dataset_subset
```

    ## An object of class Seurat 
    ## 12000 features across 385 samples within 2 assays 
    ## Active assay: RNA (2000 features, 0 variable features)
    ##  1 other assay present: peaks

``` r
# Create an hummus object from seurat object
hummus <- as(chen_dataset_subset, 'hummus_object')
```

### 1.2. Add genome and motif annotations to hummus object

Fetch genome annotations online (necessitate an internet connection).
You can also request any “EnsDB” object adapted to your data
(e.g. EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86 for human genome
annotations) or use your own genome annotations in the same format.

``` r
# get human genome annotation from EndDb data
# wrapper of Signac::GetGRangesFromEnsDb, adapting output to UCSC format
genome_annotations <- get_genome_annotations(
  ensdb_annotations = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
```

Add genome annotations to hummus/seurat object

``` r
# We add annotation to the ATAC/peak assay
# (will be used to define peak-gene links)
Signac::Annotation(hummus@assays$peaks) <- genome_annotations
rm(genome_annotations)
```

Get TF motifs from JASPAR2020 and chromVARmotifs databsases in a
motifs_db object. By default, human motifs are used. You can specify the
species you want to use with the `species` argument (e.g. species =
“mouse” for mouse). motifs_db objects contain 3 slots : \*
`motifs = "PWMatrixList"` \* `tf2motifs = "data.frame"` \*
`tfs = "NULL"` PWMatrixList is a named vector of the motif matrices,
whil tf2motifs is a correspondance table between TFs and motifs. tfs is
a named vector of the TFs. You can also use your own motifs_db object,
as long as it contains the same slots.

``` r
# Load TF motifs from JASPAR2020 and chromVARmotifs in hummus object
hummus@motifs_db <- get_tf2motifs() # by default human motifs
```

## 2. Construction of the multilayer

Since it’s a long sttep (several hours, at least to build the scATAC
layer, but also probably the scRNA layer depending of the number of
cores), you can also directly import a multilayer completed with :
`data(chen_subset_hummus)`. This object corresponds to a multilayer from
chen_dataset_subset completed. You can then go to the part 3, replacing
`hummus` by `chen_subset_hummus` in each step. \#### Compute bipartites
and add it to hummus object

### 2.1. TF - peaks bipartite reconstruction

TF - peaks bipartite is computed using the motifs_db object and the peak
assay. You can specify the assay to use to filter TFs (e.g. “RNA” if you
want to use only the TFs expressed in your dataset). If NULL, all TFs
with motifs will be used. BSGenome object is used to identify location
of motifs and intersect them with peak <br>You can also specify the name
of the bipartite that will be added to the hummus object. By default, it
will be named “tf_peak”.

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
    ##   31 TFs expressed

    ## Building motif matrix

    ## Finding motif positions

    ## Creating Motif object

    ## Warning in CreateMotifObject(data = motif.matrix, positions = motif.positions,
    ## : Non-unique motif names supplied, making unique

    ##  Adding TF info
    ##  Returning TF-peak links as bipartite object
    ## no peak layer name provided, using peak_assay name

### 2.2. Genes - peaks bipartite reconstruction

Peaks - genes bipartite is computed

``` r
hummus <- bipartite_peaks2genes(
                      hummus_object = hummus,
                      gene_assay = "RNA",
                      peak_assay = "peaks",
                      store_network = FALSE,
                      )
```

#### Compute layer networks and add it to hummus object

Each one of the three layers is computed individually.

### 2.3. Compute the TF network from OmniPath database

We currently use OmniPath R package to fetch TF interactions. You can
first specify if you want to use only the TFs expressed in your dataset
(if you have a RNA assay in your hummus object). If `gene_assay` is
NULL, all TFs with motifs will be used. <br> You can then specify which
interactions you want to keep through ‘source_target’ argument (“AND” \|
“OR”). If “AND”, only the interactions between 2 TFs that are both
present in the dataset will be kept. If “OR”, all interactions involving
at least one TF present in the dataset will be kept. <br>Finally, you
can specify the name of the multiplex and the name of the network that
will be added to the hummus object. The added network will be undirected
and unweighted since PPI and OmniPath database are not directional nor
return any weight here.

``` r
hummus <- compute_tf_network(hummus,
                            gene_assay = "RNA", # default = None ;
                                                # If a assay is provided,
                                                # only the TFs that are present
                                                # will be considered
                            verbose = 1,
                            #source_target = "OR",
                            multiplex_name = "TF",
                            tf_network_name = "TF_network")
```

    ## Computing TF network...
    ##   31 TFs expressed
    ##  TF network construction time: 1.04869 
    ## No TF-TF edges from Omnipath for the given parameters.
    ##         You can try to change the source_target parameter to 'OR' to get
    ##         TF-other protein interactions. Or try to import a network  
    ##         computed externally. Right now, a network with all TFs connected
    ##         to a fake node is created, for HuMMuS analysis.
    ##  It has no biological
    ##         meaning but will allow to run the pipeline as if no edges were present.
    ##         
    ##  Creating new multiplex :  TF

### 2.4. Compute gene network from scRNA-seq w/ GENIE3

Different methods can be used to compute the gene network. For now, only
GENIE3 is implemented in HuMMuS. You can specify which assay to use to
compute the network (`gene_assay`). <br>You can specify the number of
cores to use to compute the network. You can also specify if you want to
save the network locally (`store_network = TRUE`) or not
(`store_network = FALSE`). If you choose to save the network, you will
need to specify the output file name (`output_file`). The returned
network will be considered undirected and weighted. While GENIE3 returns
a directed network, we symmetrize it for the random walk with restart
exploration of the genes proximity.

``` r
hummus <- compute_gene_network(
              hummus,
              gene_assay = "RNA",
              method = "GENIE3",
              verbose = 1,
              number_cores = 5, # GENIE3 method can be ran
                                # parallelised on multiple cores
              store_network = FALSE, # by default : FALSE, but
                                     # each network can be saved 
                                     # when computed with hummus
              output_file = "gene_network.tsv")
```

    ## Computing gene network with  GENIE3  ...
    ##  No TFs list provided, fetching from hummus object...
    ##   31 TFs expressed
    ##  Gene network construction time: 4.797178 
    ##  Creating new multiplex :  RNA

### 2.5. Compute the peak network from scATAC-seq w/ Cicero

Different methods can be used to compute the peak network. For now, only
Cicero is implemented in HuMMuS. You can specify which assay to use to
compute the network (`peak_assay`). You can also specify the number of
cores to use to compute the network. You can also specify if you want to
save the network locally (`store_network = TRUE`) or not
(`store_network = FALSE`). If you choose to save the network, you will
need to specify the output file name (`output_file`). The returned
network will be considered undirected and weighted, since cis-regulatory
interaction and Cicero outputs are not directional.

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
    ## Mean shared cells bin-bin: 6.54129911788292
    ## Median shared cells bin-bin: 0

    ## Warning in cicero::make_cicero_cds(input_cds, reduced_coordinates =
    ## umap_coords, : On average, more than 10% of cells are shared between paired
    ## bins.

    ## [1] "Starting Cicero"
    ## [1] "Calculating distance_parameter value"

    ## Warning in estimate_distance_parameter(cds, window = window, maxit = 100, :
    ## Could not calculate sample_num distance_parameters (97 were calculated) - see
    ## documentation details

    ## [1] "Running models"
    ## [1] "Assembling connections"
    ## [1] "Successful cicero models:  4782"
    ## [1] "Other models: "
    ## 
    ## Zero or one element in range 
    ##                         8659 
    ## [1] "Models with errors:  0"
    ## [1] "Done"
    ## 
    ##  6394 peak edges with a coaccess score > 0 were found.
    ## Peak network construction time: 1.170815 Creating new multiplex :  peaks

## 3. Analyse of the multilayer and definition of GRN

### 3.1. Save the mulilayer in a classic hierarchical structure

The package used for the random walk with restart exploration
(multixrank) requires currently to save all the network files on disk.
To simplify the organisation of the file, it is possible to save
everything necessary with function `save_multilayer()`.<br> It will
create a folder (specified through `folder_name`) containing all the
files necessary to run the multixrank algorithm. The folder will contain
the following subfolders : \* **bipartite** : containing the bipartites
files \* **multiplex** : containing the multiplex sub-subfolders \*
**multiplex_1** (e.g. TF\|peak\|RNA) : containing the network file of
each layer of the multiplex \* **seed** : that will contain the seed
files (necessary to compute HuMMuS outputs later) \* **config** : that
will contain the config files (necessary to compute HuMMuS outputs
later)

``` r
save_multilayer(hummus = hummus,
                folder_name = "chen_multilayer")
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

    ## Multiplex of  1  networks with 32 features.
    ##  Networks names:  TF_network[1] "TF TF_network"

    ## Warning in dir.create(paste0(folder_name, "/", multiplex_folder, "/",
    ## multiplex_name)): 'chen_multilayer/multiplex/RNA' already exists

    ## Multiplex of  1  networks with 2000 features.
    ##  Networks names:  RNA_GENIE3[1] "RNA RNA_GENIE3"

    ## Warning in dir.create(paste0(folder_name, "/", multiplex_folder, "/",
    ## multiplex_name)): 'chen_multilayer/multiplex/peaks' already exists

    ## Multiplex of  1  networks with 4639 features.
    ##  Networks names:  peak_network_cicero[1] "peaks peak_network_cicero"

### 3.2. Retrieve target genes

With HuMMuS, inference of GRN and target gene of TFs are different
outputs. Indeed, while GRN is computed making TFs compete to regulate
genes (by random walk with restart starting from the genes and going to
the TFs), target genes are computed making genes compete to be regulated
by TFs (by random walk with restart starting from the TFs and going to
the genes).<br> For target genes output, you can specify the list of TFs
(`tf_list`) to use as seed (if NULL by default, all TFs will be used as
seed). Only the links between the seed TFs and the genes will be
computed. You can also specify the list of genes to use. Only the score
of the genes present in the network and the `gene_list` will be
returned.

``` r
target_genes <- define_target_genes(
  hummus,
  multilayer_f = "chen_multilayer",
  njobs = 5
  )
```

``` r
head(target_genes)
```

    ## [1] layer      path_layer score      tf         gene      
    ## <0 rows> (or 0-length row.names)

### 3.3. Define GRN

The GRN is defined using the multixrank algorithm. It requires to have
the hummuspy python package installed (pip install hummuspy).<br> The
function `define_grn()` will call hummuspy and compute random walk with
restart from each gene to find the most probable TFs regulating it. It
will return a dataframe containing the TF-gene links with their
associated probability.<br> This can be parallelised using the njobs
argument. You can also specify the list of genes to use as seed (if NULL
by default, all genes will be used as seed). Only the links between the
seed genes and the TFs will be computed.<br> You can also specify the
list of TFs to use. If not NULL, the score of the TFs present in the
network and the `tf_list` will be returned. <br> In addition, you can
also specify the location of each file if you didn’t use the default
option when saving the multilayer. Finally, you can also specify if you
want to return or to save the output. If you choose to save the output,
you will need to specify the output file name. Output returned is a
dataframe containing TF-gene links with their associated probability.

**Similar options are available for the other outputs (enhancers,
binding regions and target genes).**

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

    ## [1] layer      path_layer score      gene       tf        
    ## <0 rows> (or 0-length row.names)

### 3.4. Retrieve enhancers

For enhancers output, you can specify the list of genes (`gene_list`) to
use as seed (if NULL by default, all genes will be used as seed). Only
the links between the seed genes and the peaks will be computed. <br>
You can also specify the list of peaks to use. Only the score of the
peaks present in the network and the `peak_list` will be returned.

``` r
enhancers <- define_enhancers(
  hummus,
  multilayer_f = "chen_multilayer",
  njobs = 5
  )
```

``` r
head(enhancers)
```

    ## [1] layer      path_layer score      gene       peak      
    ## <0 rows> (or 0-length row.names)

### 3.5. Retrieve binding regions

For binding regions output, you can specify the list of TFs (`tf_list`)
to use as seed (if NULL by default, all TFs will be used as seed). Only
the links between the seed TFs and the peaks will be computed. You can
also specify the list of peaks to use. Only the score of the peaks
present in the network and the `peak_list` will be returned.

``` r
binding_regions <- define_binding_regions(
  hummus,
  multilayer_f = "chen_multilayer",
  njobs = 5
  )
```

``` r
head(binding_regions)
```

    ## [1] layer      path_layer score      tf         peak      
    ## <0 rows> (or 0-length row.names)
