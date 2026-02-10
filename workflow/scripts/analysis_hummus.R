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

# Creating the reports folder
rlog::log_info("[Step 3 | analysis_hummus] - Creating the reports folder")
reports_folder = file.path(wrkdir, "reports")
dir.create( reports_folder, showWarnings = FALSE )

# -------------------------------------------------------------
# Loading saved hummus object
rlog::log_info("[Step 3 | analysis_hummus] - Loading saved hummus object")
obj_path = file.path(wrkdir, "hummus_object.rds")
hummus = readRDS(obj_path)

cores = 4
if( ! is.null(config$resources$cores) ){
  cores = config$resources$cores
}

# -------------------------------------------------------------
# Performing the analysis stated in the configuration file
rlog::log_info("[Step 3 | analysis_hummus] - Performing the analysis stated in the configuration file")
cnf_analysis = config$analysis
if( ! is.null(cnf_analysis) ){

  if( ! is.null(cnf_analysis$grn) ){
    seed_gene = cnf_analysis$grn$seeds$gene
    seed_tf = cnf_analysis$grn$seeds$tf
    fname = cnf_analysis$grn$filename
    if( is.null(fname) ){
      fname = "ranked_grn.tsv"
    }
    fname = file.path(reports_folder, basename(fname) )
    rlog::log_info("[Step 3 | analysis_hummus] - Performing analysis >> grn")
    out <- define_grn( hummus, gene_list = seed_gene, tf_list = seed_tf, save=TRUE, output_f = fname, njobs = cores )
  }

  if( ! is.null(cnf_analysis$target_gene) ){
    seed_gene = cnf_analysis$target_gene$seeds$gene
    seed_tf = cnf_analysis$target_gene$seeds$tf
    fname = cnf_analysis$target_gene$filename
    if( is.null(fname) ){
      fname = "ranked_target_gene.tsv"
    }
    fname = file.path(reports_folder, basename(fname) )
    rlog::log_info("[Step 3 | analysis_hummus] - Performing analysis >> target genes")
    out <- define_target_genes( hummus, gene_list = seed_gene, tf_list = seed_tf, save=TRUE, output_f = fname, njobs = cores )
  }

  if( ! is.null(cnf_analysis$enhancer) ){
    seed_gene = cnf_analysis$enhancer$seeds$gene
    seed_tf = cnf_analysis$enhancer$seeds$tf
    fname = cnf_analysis$enhancer$filename
    if( is.null(fname) ){
      fname = "ranked_enhancer.tsv"
    }
    fname = file.path(reports_folder, basename(fname) )
    rlog::log_info("[Step 3 | analysis_hummus] - Performing analysis >> enhancers")
    out <- define_enhancers( hummus, gene_list = seed_gene, tf_list = seed_tf, save=TRUE, output_f = fname, njobs = cores )
  }

  if( ! is.null(cnf_analysis$binding_region) ){
    seed_gene = cnf_analysis$binding_region$seeds$gene
    seed_tf = cnf_analysis$binding_region$seeds$tf
    fname = cnf_analysis$binding_region$filename
    if( is.null(fname) ){
      fname = "ranked_binding_region.tsv"
    }
    fname = file.path(reports_folder, basename(fname) )
    rlog::log_info("[Step 3 | analysis_hummus] - Performing analysis >> binding region")
    out <- define_binding_regions( hummus, gene_list = seed_gene, tf_list = seed_tf, save=TRUE, output_f = fname, njobs = cores )
  }
}