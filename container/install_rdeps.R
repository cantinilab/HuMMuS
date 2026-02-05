chooseCRANmirror(ind = 1)

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

install.packages("/opt/HuMMuS", repos = NULL, type = "source")
