
## Set-up environment

#setwd('../Desktop/PhD/Multilayer networks/example_HMLN/')

# GENIE3 to infer gene links
NUM_CORES <- 50       # Number of cores for network construction with GENIE3

# cicero to infer peak links
NUMBER_CLUSTER <- 50           # Default 50
NUMBER_CLUSTER_low <- 10       # We use this instead of 50 for hESC_Liu since it has onlly 72 cells
NUMBER_SAMPLE <- 100           # Default 100 # never used!
MY_SEED <- 2021

# Parameters for determining peak to gene links
PEAK_TO_GENE_METHOD = 'Signac'
UPSTREAM = 500                 # Distance upstream of the gene to consider as potential regulatory region
DOWNSTREAM = 500               # Distance downstream of the gene to consider as potential regulatory region
ONLY_TSS = TRUE                # Measure distance from the TSS (TRUE) or from the entire gene body (FALSE)


## Dependencies

# Bipartites: TF-peak and peak-gene links
library(Seurat)
library(Signac)
#devtools::install_github('quadbiolab/Pando')
#library(Pando)
library(JASPAR2020)
library(chromVARmotifs)
library(plyr)
library(dplyr)
library(tibble)
library(tidyr)
library(data.table)
library(TFBSTools)
library(stringr)
library(GENIE3)

# Mouse
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(BSgenome.Mmusculus.UCSC.mm10)
#BiocManager::install("EnsDb.Mmusculus.v79")
#library(BSgenome.Mmusculus.UCSC.mm10)
#library(EnsDb.Mmusculus.v79)

# Human
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#BiocManager::install("EnsDb.Hsapiens.v86")
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)

# GENIE3 to infer gene links
#BiocManager::install("GENIE3")
library(GENIE3)
#library(doParallel)
#library(doRNG)

# cicero to infer peak links
library(reshape2)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(c("Gviz", "GenomicRanges", "rtracklayer")) ############### not sure if needed???? to test if it also runs without
#devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3") # to have the good version of Cicero
library(cicero)
library(monocle3)
#library(dplyr) # installed above
#library(tidyr) # installed above

# remap
#library(reshape2) # installed above
#install.packages("Hmisc")
#library(Hmisc)

# Download TF - TF links from OmnipathR database
library(OmnipathR)

# evaluation visualization
#library(ggplot2)
#library(RColorBrewer)
#library(gridExtra)

