#!/usr/bin/env python
import argparse
import os
import mudata as md
import perturbo
import scvi

scvi.settings.seed = 0

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Run PerTurbo on a single dataset")
parser.add_argument("-f", "--file", required=True, help="Path to a MuData object")
parser.add_argument("-o", "--output", required=True, help="Name of output .csv file")
parser.add_argument("-m", "--min_umi_count", default=1000, help="Minimum UMI count")
parser.add_argument("-r", "--rna_modality", default="rna", help="RNA modality name")
parser.add_argument(
    "-g", "--guide_modality", default="grna", help="guide RNA (gRNA) modality name"
)
parser.add_argument("-bk", "--batch_key", help="Batch ID key")
parser.add_argument("-l", "--library_size_key", help="Library size key")
parser.add_argument("-e", "--epochs", default=50, help="Number of epochs")

parser.add_argument(
    "-b", "--batch_size", default=512, help="Training batch size (# of cells)"
)
args = parser.parse_args()

if os.path.splitext(args.file)[1] != ".h5mu":
    raise Exception("Error: file must have .h5mu extension")

batch_size = args.batch_size

# Load and register data
mdata = md.read(args.file)
perturbo.PERTURBO.setup_mudata(
    mdata,
    batch_key=args.batch_key,
    library_size_key=args.library_size_key,
    continuous_covariates_keys=["percent_mito"],
    guide_by_element_key="element_targeted",
    modalities={
        "rna_layer": args.rna_modality,
        "perturbation_layer": args.guide_modality,
    },
)

# Set up and train model
model = perturbo.PERTURBO(mdata, likelihood="nb")
model.view_anndata_setup()
model.train(args.epochs, lr=0.01, batch_size=batch_size)

element_effects_df = model.get_element_effects()
mdata.uns["test_results"] = element_effects_df
mdata.write_h5mu(args.output)
