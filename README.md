![Build](https://github.com/cantinilab/HuMMuS/____/badge.svg?branch=main)

# HuMMuS <img src="Figures/hummus_logo.png" align="right" width="180"/>
## Heterogeneous Multilayer network for Multi-omics Single-cell data 

HuMMuS exploits multi-omics single-cell measurements to infer numerous regulatory relationships.
Beside classical Gene Regulatory Networks (GRN), HuMMus proposes enhancer prediction, binding regions prediction and target genes of specific transcription factors.

#### **scRNA + scATAC** 
Like most of the current state-of-the-art method to infer Gene Regulatory Networks, we popose a minimal version of HuMMuS based on scRNA-seq + scATAC-seq data (paired or **unpaired**).

#### **Use of additional modalities**
HuMMuS has been developed to be extendable to any additional biological modality of interest.
It is then possible to add any additional network to an already existing modality (forming multiplex where the same nodes are connected through different networks (e.g. prior-knowledge network and data-driven network)) or from a new one (e.g. adding epigenetic or proteomic networks).


## Installation

HuMMuS is for now ready only in R but requires some python dependencies (hummuspy).

#### HuMMuS python depency
Python package hummuspy should preferably be installed using pip (from the terminal in a conda environment for e.g)

```r
conda create -n hummuspy_env
conda activate hummuspy_env
pip install hummuspy
```

Alternatively, you can also install it directly from R using the reticulate package:
```r
library(reticulate)
py_install("hummuspy", envname = "r-reticulate", method="auto")
```

Before running HuMMuS, if you're using multiple conda environment you need to make sure to that reticulate points toward the one where hummuspy is installed. You can precise it at the beginning of your code :

```r
library(reticulate)
# Using a specific conda environment
envname = "hummuspy_env" # or "r-reticulate" for e.g.
use_condaenv(envname, required = TRUE)
```

### HuMMuS R package
Core R package can be installed directly from R:
```r
devtools::install_github("cantinilab/HuMMuS") 
```

For more details on hom to setup the reticulate connection,
see: https://rstudio.github.io/reticulate

## Tutorials/Vignettes

The tutorial and vignettes proposed will be listed here. For now, we propose a vignette to illustrate the most standard use of HuMMuS.
* **Infer a gene regulatory network and other outputs from unpaired/paired scRNA+scATAC data**: This [vignette](https://github.com/cantinilab/HuMMuS/blob/main/examples/chen_grn.md) 
shows the application of HuMMuS to the Chen dataset, used in the benchmark of HuMMuS publication [paper](__preprint_links__).

## Data accessibility

To reproduce HuMMuS results presented in the manuscript, preprocessed data are accessible [here](https://figshare.com/account/home#/projects/168899)
<br> For quick tests, the Chen dataset preprocessed is accessible directly through the package as a Seurat object: `load(chen_data)`.

## Cite us
Trimbour R., Deutschmann I. M., Cantini L. Molecular mechanisms reconstruction from single-cell multi-omics data with HuMMuS. bioXriv (2023). doi: 







