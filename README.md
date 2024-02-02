# HuMMuS <img src="man/figures/hummus_logo.png" align="right" width="180"/>
### Heterogeneous Multilayer network for Multi-omics Single-cell data

HuMMuS exploits multi-omics single-cell measurements to infer numerous regulatory mechanisms.
Inter-omics (e.g. peak-gene, TF-peak) and intra-omics interactions (e.g. peak-peak, gene-gene, TF-TF) are considered to capture both regulatory interactions and macromolecule cooperations.

## Overview

The current outputs available from HuMMuS are

* gene regulatory networks (GRNs)
* enhancers
* TF - DNA binding regions
* TF - target genes.

#### [Read our preprint](https://www.biorxiv.org/content/10.1101/2023.06.09.543828v1) for more details !
<img src="man/figures/Fig_0001.jpg" align="center" width="1000"/>

### **scRNA + scATAC**
Like most of the current state-of-the-art methods to infer GRN, we propose a minimal version of HuMMuS based on scRNA-seq + scATAC-seq data (paired or **unpaired**).

### **Use of additional modalities**
HuMMuS has been developed to be extendable to any additional biological modality of interest.
It is then possible to add any additional network to an already **existing modality** (e.g. both prior-knowledge network and data-driven network of genes), or from a **new modality** (e.g. adding epigenetic or proteomic networks).
<br>_For now, such personalisation requires to use directly some hummuspy (python package) functions at the end of the pipeline and write some configuration files manually. It will be simplified soon !_

## Tutorials/Vignettes

* [**Infer a gene regulatory network and other outputs from unpaired/paired scRNA+scATAC data**](https://github.com/cantinilab/HuMMuS/blob/main/vignettes/chen_grn.md) shows the application of HuMMuS to the Chen dataset, used in the benchmark of [HuMMuS publication](https://www.biorxiv.org/content/10.1101/2023.06.09.543828v1).

## Installation
HuMMuS is for now ready only in R but requires some python dependencies (hummuspy).

### HuMMuS python dependency
Python package **hummuspy** should preferably be installed using pip (from the terminal in a conda environment for e.g)
```r
conda create -n hummuspy_env python
conda activate hummuspy_env
pip install hummuspy
```

Alternatively, you can also install it directly from R using the reticulate package:
```r
library(reticulate)
py_install("hummuspy", envname = "r-reticulate", method="auto")
```

### HuMMuS R package
Core R package can be installed directly from R:
```r
devtools::install_github("cantinilab/HuMMuS")
```

Before running HuMMuS, if you're using multiple conda environment you need to make sure to that reticulate points toward the one where hummuspy is installed. You can precise it at the beginning of your code :
```r
library(reticulate)
# Using a specific conda environment
envname = "hummuspy_env" # or "r-reticulate" for e.g.
use_condaenv(envname, required = TRUE)
```
For more details on how to setup the reticulate connection,
see: https://rstudio.github.io/reticulate

### scATAC processing
To compute directly the scATAC data with HuMMuS, we currently only propose to use [Cicero](https://cole-trapnell-lab.github.io/cicero-release/docs_m3/). It requires the version running with [Monocle3](https://cole-trapnell-lab.github.io/monocle3/).
You then need to install both [Monocle3](https://cole-trapnell-lab.github.io/monocle3/docs/installation/), and Cicero:

```r
devtools::install_github("cole-trapnell-lab/monocle3")
devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")
```
*If you encounter some troubles with Monocle3 installation, on ubuntu you can try to run: `sudo apt-get install libgdal-dev libgeos-dev libproj-dev`. You can also go on [their github page](https://github.com/cole-trapnell-lab/monocle3/issues) for more help. Having a previous version of Monocle (1 or 2) still in your R session can cause some troubles. If you encounter some even after restarting your R session,, try to `remove.packages("monocle")` before to reinstall both Monocle**3** and Cicero*

## Data accessibility

To reproduce HuMMuS results presented in the manuscript, preprocessed data [are accessible here](https://figshare.com/account/home#/projects/168899)
<br> For quick tests, the Chen dataset preprocessed is accessible directly through the package as a Seurat object: `load(chen_data)`.

## Cite us
Trimbour R., Deutschmann I. M., Cantini L. Molecular mechanisms reconstruction from single-cell multi-omics data with HuMMuS. bioXriv (2023). doi: https://doi.org/10.1101/2023.06.09.543828
