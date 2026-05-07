#!/usr/bin/env Rscript
library(reticulate)
reticulate::use_condaenv('rhummus_env')
hummuspy <- reticulate::import("hummuspy")
library("optparse")
library("jsonlite")
library("HuMMuS")
 
option_list = list(
  make_option(c("-d", "--work_dir"), type="character", default=NULL, 
              help="Working directory", metavar="character"),
  make_option(c("-f", "--configuration_file"), type="character", default=NULL, 
              help="Configuration file", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
config = fromJSON(opt$configuration_file)
wrkdir = opt$work_dir

# -------------------------------------------------------------
# Creating the downloads folder
rlog::log_info("[Step 2 | generate_networks_layers] - Loading Seurat assays data")
downloads_folder = file.path(wrkdir, "downloads")
dir.create( downloads_folder, showWarnings = FALSE )

# Loading saved hummus object
rlog::log_info("[Step 2 | generate_networks_layers] - Loading saved hummus object")
obj_path = file.path(wrkdir, "hummus_object.rds")
hummus = readRDS(obj_path)

# -------------------------------------------------------------
# Computing layers
rlog::log_info("[Step 2 | generate_networks_layers] - Computing layers")

rlog::log_info("[Step 2 | generate_networks_layers] - Computing layers >> TF-motifs")
hummus@motifs_db <- get_tf2motifs( download_folder = downloads_folder )

rlog::log_info("[Step 2 | generate_networks_layers] - Computing layers >> TF-peaks")
hummus <- bipartite_tfs2peaks( hummus_object = hummus, tf_expr_assay = "RNA", peak_assay = "peaks", tf_multiplex_name = "TF", genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38 )

rlog::log_info("[Step 2 | generate_networks_layers] - Computing layers >> peaks-genes")
hummus <- bipartite_peaks2genes( hummus_object = hummus, gene_assay = "RNA", peak_assay = "peaks", store_network = FALSE )

# -------------------------------------------------------------
# Generating multiplex networks
rlog::log_info("[Step 2 | generate_networks_layers] - Generating multiplex networks")

rlog::log_info("[Step 2 | generate_networks_layers] - Generating multiplex networks >> TF")
method_tf = config$network_methods$tf
if( is.null(method_tf) ){
  method_tf = "Omnipath"
}
hummus <- compute_tf_network(hummus, gene_assay = "RNA", method = method_tf, verbose = 1, multiplex_name = "TF", tf_network_name = "TF_network")

rlog::log_info("[Step 2 | generate_networks_layers] - Generating multiplex networks >> gene")
method_gene = config$network_methods$gene
if( is.null(method_gene) ){
  method_gene = "GENIE3"
}
hummus <- compute_gene_network( hummus, gene_assay = "RNA",  method = method_gene,  verbose = 1, number_cores = 5)

rlog::log_info("[Step 2 | generate_networks_layers] - Generating multiplex networks >> peak")
method_peak = config$network_methods$peak
if( is.null(method_peak) ){
  method_peak = "cicero"
}
hummus <- compute_atac_peak_network(hummus, atac_assay = "peaks", method = method_peak, verbose = 1, genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38) 

# -------------------------------------------------------------
# Storing and updating hummus object
store_update_hummus_object_wrapper(hummus)

# Storing multilayer network
rlog::log_info("[Step 2 | generate_networks_layers] - Storing multilayer network")
save_multilayer(hummus)