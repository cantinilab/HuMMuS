![Build](https://github.com/cantinilab/HuMMuS/____/badge.svg?branch=main)

# HuMMuS <img src="Figures/hummus_logo.png" align="right" width="180"/>
### Heterogeneous Multilayer network for Multi-omics Single-cell data 

HuMMuS exploits multi-omics single-cell measurements to infer numerous regulatory relationships.
Beside classical Gene Regulatory Networks (GRN), HuMMus propose enhancer predictions, binding regions prediction and even broad communities detection.

#### **scRNA + scATAC** 
Like most of the current state-of-the-art method to infer Gene Regulatory Networks, we popose a minimal version of HuMMuS based on scRNA-seq + scATAC-seq data (paired or **unpaired**).

#### **Use of additional modalities**
HuMMuS has been developed to be extendable to any additional biological modality of interest.
It is then possible to add any additional network to an already existing modality (forming multiplex where the same nodes are connected through different networks (e.g. prior-knowledge network and data-driven network)) or from a new one (e.g. adding epigenetic or proteomic networks).


## Installation

HuMMuS is for now ready only in R but it requires some python dependencies (hummuspy).

### Python depencency : hummuspy
Python package hummuspy should preferably be installed using pip (from the terminal in your conda environment for e.g)
```r
pip install hummuspy

Alternatively, you can also install it directly from R using the reticulate package:
library(reticulate)
envname = "hummus_env"
py_install("hummuspy", envname = envname, method="auto")
```

### HuMMuS R package
Core R package can be installed directly from R:
```r
devtools::install_github("cantinilab/hummus") 
```

### Reticulate configuration
Before running HuMMuS, you need to make sure 'reticulate' is using the correct python binary (or conda environment). You can precise it at the beginning of your code :
```r
library(reticulate)

# Using a specific python binary
use_python("/home/user/python", required=TRUE)

# Using a specific conda environemnt
envname = "hummus_env"
use_condaenv(envname, required = TRUE)
```

For more details on hom to setup the reticulate connection,
see: https://rstudio.github.io/reticulate


## Tutorials/Vignettes

The tutorial and vignettes proposed will be listed here. For now, we propose a vignette to illustrate the most standard use of HuMMuS.
* **Infer a gene regulatory network from unpaired/paired scRNA+scATAC data** The data correspond to a subset of the Chen dataset, which was part of the benchmark analysed in the [paper] __preprint_links__




