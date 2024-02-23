## Helper script to create simulation muData file from inference example

library(MuData)
library(dplyr)
library(tidyr)

## Replace with NextFlow input and output
infile <- "/mnt/shared/inference/gasperini_inference_input.h5mu"
outfile <- "output/gasperini_simulation_input.h5mu"

# load muData file
mu <- readH5MU(infile)

# pairs to perturb
pairs_to_perturb <- data.frame(gene_id = c("ENSG00000136856", "ENSG00000136925"),
                      intended_target_name = c("candidate_enh_3", "candidate_enh_2"),
                      effect_size = 0.5)

# add pairs_to_test back to muData object
mu@metadata$pairs_to_test <- as.list(pairs_to_perturb)

# save muData object to output file
writeH5MU(mu, file = outfile)
