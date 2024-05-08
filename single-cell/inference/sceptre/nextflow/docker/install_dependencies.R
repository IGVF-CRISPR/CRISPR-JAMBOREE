# The purpose of this script is to create a renv.lock file that contains all the
# R packages on which sceptreIGVF and their versions. It should be rerun only if
# one wishes to update the R package versions in renv.lock. Generally, we want
# to freeze the R package versions in renv.lock for reproducibility purposes, so
# this script will seldom be rerun.

install.packages("BiocManager", repos = "http://cran.us.r-project.org")
BiocManager::install(version = "3.19")
BiocManager::install("MultiAssayExperiment")
BiocManager::install("SingleCellExperiment")
BiocManager::install("rhdf5")

install.packages("remotes")
remotes::install_github("ilia-kats/MuData")
remotes::install_github("katsevich-lab/sceptre")
remotes::install_github("IGVF-CRISPR/sceptreIGVF")

install.packages("renv")
renv::settings$ignored.packages(c("remotes"), persist = FALSE)
renv::settings$bioconductor.version("3.19")
renv::init()
