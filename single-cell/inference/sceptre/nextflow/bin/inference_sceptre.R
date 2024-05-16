#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# obtain the command line arguments
mudata_fp <- args[1]
side <- args[2]
grna_integration_strategy <- args[3]
resampling_approximation <- args[4]
control_group <- args[5]
resampling_mechanism <- args[6]
formula_object <- args[7]

# process formula object
if (!identical(formula_object, "default")) {
  formula_object <- stats::formula(formula_object)
}

# read MuData
mudata_in <- MuData::readH5MU(mudata_fp)

# run sceptre inference
mudata_out <- sceptreIGVF::inference_sceptre(
  mudata = mudata_in,
  side = side,
  grna_integration_strategy = grna_integration_strategy,
  resampling_approximation = resampling_approximation,
  control_group = control_group,
  resampling_mechanism = resampling_mechanism,
  formula_object = formula_object
)

# write MuData
MuData::writeH5MU(object = mudata_out, file = 'mudata_out.h5mu')