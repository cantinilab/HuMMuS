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

reticulate::use_condaenv('rhummus_env')

# Loading Seurat assays data
rlog::log_info("[Step 1 | process_data] - Loading Seurat assays data")
assay_data = load_first_object( config$assays_rda_path )

# -------------------------------------------------------------
# Initialize Hummus object
rlog::log_info("[Step 1 | process_data] - Initializing Hummus object")
hummus <- Initiate_Hummus_Object(assay_data, multilayer_folder = wrkdir )

# -------------------------------------------------------------
# Annotating peaks
rlog::log_info("[Step 1 | process_data] - Annotating peaks")
genome_annotations <- get_genome_annotations(  ensdb_annotations = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86 )
Signac::Annotation(hummus@assays$peaks) <- genome_annotations
rm(genome_annotations)

# -------------------------------------------------------------
# Saving object
rlog::log_info("[Step 1 | process_data] - Saving hummus object")
store_update_hummus_object_wrapper(hummus)