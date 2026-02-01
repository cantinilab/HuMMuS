conda create -n hummus_env python=3.10 joblib pandas numpy tqdm pyyaml matplotlib scipy distributed dask bokeh rich pyarrow r-base
conda activate hummus_env
conda install bioconda::multixrank
conda install r::r-reticulate
conda install r::r-seurat
conda install bioconda::r-signac
conda install -c conda-forge r-biocmanager

Call Rscript install_deps.R , this file contains:

BiocManager::install(update=TRUE, ask=FALSE)
BiocManager::install("sparseMatrixStats")
BiocManager::install("TFBSTools")



------------------
conda create -n hummus_env python=3.10 pandas numpy r-base
conda activate hummus_env
conda install r-devtools

sudo apt install tar
sudo ln -s /usr/bin/tar /bin/gtar

Call Rscript install_deps.R , this file contains:
library("devtools")                                                                                                             
devtools::install_github("cantinilab/HuMMuS", ref="dev_SeuratV5")
