args <- commandArgs(trailingOnly = TRUE)
mudata_fp <- args[1]
output_fp <- args[2]

# 1. library the packages
library(MuData)
library(sceptreIGVF)

# 2. load the mudata
mudata <- readH5MU(file = mudata_fp)
# data("mudata_guide_assignment")
# mudata <- mudata_guide_assignment
# output_fp <- "example.h5mu"

# 3. call sceptre grna assignment on the mu data
mudata_updated <- assign_grnas_sceptre(mudata = mudata)

# 4. write mudata
writeH5MU(mudata_updated, file = output_fp)
