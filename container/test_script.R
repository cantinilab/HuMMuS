/*
Environment

sudo apt-get install libudunits2-dev libgdal-dev libgeos-dev libproj-dev 
sudo apt-get install zlib1g-dev libabsl-dev

conda create --name rhummus_env python=3.10 r-base=4.4.3 -y
conda install h5py
conda install r-reticulate r-jsonlite r-devtools r-signac r-grr r-doparallel r-dorng r-monocle3

conda install bioconda::bioconductor-tfbstools bioconda::bioconductor-chromvar bioconda::bioconductor-motifmatchr bioconda::bioconductor-singlecellexperiment bioconda::bioconductor-biovizbase 
conda install omnipath bioconda::bioconductor-singlecellexperiment
conda install numpy=1.26.4

chooseCRANmirror(ind = 1)

Install the rest in r cmd
library(reticulate)
reticulate::use_condaenv('rhummus_env')

library(devtools)
devtools::install_github("cantinilab/HuMMuS", ref="dev_SeuratV5", upgrade = "never")
BiocManager::install("OmnipathR")

options(timeout = 1000)
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

remove.packages("cicero")
#devtools::install_github("cole-trapnell-lab/monocle3", upgrade = "never")
devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3", upgrade = "always")
devtools::install_github("coletl/coler")

*/


folder = "/mnt/yasdata/home/yasmmin/Dropbox/portfolio_2025/hummus_cantinilab_issues/HuMMuS"
src_folder = file.path( folder, "R" )
sample_data = file.path( folder, "data" )

files = list.files( src_folder )
for (f in files) {
    source( file.path(src_folder, f) )
}
------------------
# Testing save_config
library(reticulate)
reticulate::use_condaenv('rhummus_env')
hummuspy <- reticulate::import("hummuspy")
#in playground folder

folder = "/mnt/yasdata/home/yasmmin/Dropbox/portfolio_2025/hummus_cantinilab_issues/HuMMuS"
sample_data = file.path( folder, "data" )
load(file = file.path(sample_data, "chen_dataset_subset.rda") )

src_folder = file.path( folder, "R" )
files = list.files( src_folder )
for (f in files) {
    source( file.path(src_folder, f) )
}

hummus <- Initiate_Hummus_Object(chen_dataset_subset, multilayer_folder="test_wrapper_init_dir")
grn <- define_grn( hummus, multilayer_folder = "../chen_multilayer", njobs = 5 )

# Flow ---------------

# rule 1 - prepare_object
load(file='./HuMMuS/data/chen_dataset_subset.rda')
hummus <- Initiate_Hummus_Object(chen_dataset_subset, multilayer_folder = "hmexec_dir")

genome_annotations <- get_genome_annotations(  ensdb_annotations = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86 )
Signac::Annotation(hummus@assays$peaks) <- genome_annotations
rm(genome_annotations)

# rule 2 - compute networks and layers
hummus@motifs_db <- get_tf2motifs( download_folder="./" )
hummus <- bipartite_tfs2peaks( hummus_object = hummus, tf_expr_assay = "RNA", peak_assay = "peaks", tf_multiplex_name = "TF", genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38 )
hummus <- bipartite_peaks2genes( hummus_object = hummus, gene_assay = "RNA", peak_assay = "peaks", store_network = FALSE )

hummus <- compute_tf_network(hummus, gene_assay = "RNA", method = "Omnipath", verbose = 1, multiplex_name = "TF", tf_network_name = "TF_network")
hummus <- compute_gene_network( hummus, gene_assay = "RNA",  method = "GENIE3",  verbose = 1, number_cores = 5,  store_network = TRUE, output_file = "gene_network.tsv")
hummus <- compute_atac_peak_network(hummus, atac_assay = "peaks", verbose = 1, genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,  store_network = TRUE, output_file = "peak_network.tsv") 

save_multilayer(hummus, folder_name = "chen_multilayer")

# rule 3 - analysis
ATF2_genes <- define_target_genes(  hummus, tf_list = list("ATF2"), multilayer_f = "chen_multilayer", njobs = 1)
  
target_genes <- define_target_genes( hummus, multilayer_folder = "chen_multilayer", njobs = 1)
grn <- define_grn( hummus, multilayer_folder = "chen_multilayer", njobs = 5 )
enhancers <- define_enhancers( hummus, gene_list = list("ATF2"),  multilayer_folder = "chen_multilayer", njobs = 1 )
binding_regions <- define_binding_regions( hummus, multilayer_folder = "chen_multilayer", njobs = 1 )

# -------- items snakemake parameter file
njobs = 4
working_dir = 'wrkdir_hummus'
assays_rda_path = 'chen_dataset_subset.rda'
network_methods = { 'tf': 'Omnipath', 'gene': 'GENIE3', 'peak': 'cicero' }

conda install rpy2 pyreadr - do not work
hm = readRDS('./playground/test_wrapper_init_dir/hummus_object.rds')

working_dir: "wrkdir_hummus"
assays_rda_path = 'chen_dataset_subset.rda'

resources:
    cores: 4
    
assays_rda_path: 'chen_dataset_subset.rda'
network_methods:
    tf: "Omnipath"
    gene: "GENIE3"
    peak: "cicero"
    
analysis:
    grn:
        seeds:
            gene: "all"
            tf: "all"
        save_flag: true
        name: 'ranked_grn.tsv'
    target_gene:
        seeds:
            gene: "all"
            tf: "all"
        save_flag: true
        name: 'ranked_gene.tsv'
    enhancer:
        seeds:
            gene: "all"
            tf: "all"
        save_flag: true
        name: 'ranked_enhancer.tsv'
    binding_region:
        seeds:
            gene: "all"
            tf: "all"
        save_flag: true
        name: 'ranked_binding_region.tsv'

# ----------------------------- playground ---------------------------------
# R version does not create config yml file

multiplex_names = NULL
  bipartites_names = NULL
  config_name = "config.yml"
  config_folder = "config"
  tf_multiplex = "TF"
  atac_multiplex = "peaks"
  rna_multiplex = "RNA"
  multilayer_folder = "chen_multilayer"
  gene_list = NULL
  tf_list = NULL
  save = FALSE
  output_f = NULL
  return_df = TRUE
  suffix_bipartites = ".tsv"
  njobs = 1
  output_type = 'grn'
  hummus_object = hummus
  
  # Format multiplexes names
  multiplexes_dictionary <- format_multiplex_names(
    hummus_object,
    multiplex_names = multiplex_names)
  # Format bipartites names
  bipartites_dictionary <- format_bipartites_names(
    hummus_object,
    bipartites_names = bipartites_names,
    suffix_bipartites = suffix_bipartites)

  source_python("save_config_testa.py")

# ro solve the problem "Error: Unable to access object (object is from previous session and is now invalid)" , delete the function objects imported from python file:
remove(test)
remove(get_output_from_dictsv3)
remove(define_grn_from_configv2)

  # define target_genes with hummuspy function
  output <- get_output_from_dictsv3(
    output_request = output_type,
    multilayer_folder = multilayer_folder,
    multiplexes_list = multiplexes_dictionary,
    bipartites_list = bipartites_dictionary,
    gene_list = gene_list,
    tf_list = tf_list,
    config_filename = config_name,
    config_folder = config_folder,
    output_f = output_f,
    tf_multiplex = tf_multiplex,
    peak_multiplex = atac_multiplex,
    rna_multiplex = rna_multiplex,
    update_config = TRUE,
    save = save,
    return_df = return_df,
    njobs = njobs,
    save_configfile=TRUE)

# in playground folder
library(reticulate)
reticulate::use_condaenv('rhummus_env')
hummuspy <- reticulate::import("hummuspy")
load('../.RData')


output <- hummuspy$core_grn$get_output_from_dicts(
    output_request = output_type,
    multilayer_folder = multilayer_folder,
    multiplexes_list = multiplexes_dictionary,
    bipartites_list = bipartites_dictionary,
    gene_list = gene_list,
    tf_list = tf_list,
    config_filename = config_name,
    config_folder = config_folder,
    output_f = output_f,
    tf_multiplex = tf_multiplex,
    peak_multiplex = atac_multiplex,
    rna_multiplex = rna_multiplex,
    update_config = TRUE,
    save = save,
    return_df = return_df,
    njobs = njobs,
    save_configfile=TRUE)
    
    
