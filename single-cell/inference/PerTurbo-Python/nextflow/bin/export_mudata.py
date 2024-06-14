#!/usr/bin/env python
# Exports matrices and dataframes represented in a MuData object to .mtx and .csv files

import mudata as md
import sys
from scipy.io import mmwrite
from scipy.sparse import issparse
import pandas as pd

mdata = md.read(sys.argv[1])
mmwrite("rna.mtx", mdata["rna"].X, field="integer")
mmwrite("grna.mtx", mdata["grna"].X, field="integer")

mdata["rna"].var.index.name = "gene_id"
mdata["grna"].var.index.name = "guide_id"
mdata["rna"].obs.index.name = "cell_id"

mdata["rna"].var.to_csv("gene_ids.csv")
mdata["grna"].var.to_csv("guide_ids.csv")
mdata["rna"].obs.to_csv("cell_covariates.csv")


if "element_targeted" in mdata["grna"].varm:
    guide_by_element = mdata["grna"].varm["element_targeted"]
    if isinstance(guide_by_element, pd.DataFrame):
        guide_element_df = guide_by_element
    else:
        if issparse(guide_by_element):
            guide_by_element = guide_by_element.todense()
        guide_element_df = pd.DataFrame(
            guide_by_element,
            columns=mdata["grna"].uns["elements"],
            index=mdata["grna"].var_names,
        )
    guide_element_pairs = guide_element_df.stack().reset_index()
    guide_element_pairs.columns = ["guide_id", "element", "value"]
    guide_element_pairs.query("value>0").to_csv("guide_element_pairs.csv", index=False)


if "element_tested" in mdata["rna"].varm:
    gene_by_element = mdata["rna"].varm["element_tested"]
    if isinstance(gene_by_element, pd.DataFrame):
        gene_element_df = gene_by_element
    else:
        if issparse(gene_by_element):
            gene_by_element = gene_by_element.todense()
        gene_element_df = pd.DataFrame(
            gene_by_element,
            columns=mdata["rna"].uns["elements"],
            index=mdata["rna"].var_names,
        )
    gene_element_pairs = gene_element_df.stack().reset_index()
    gene_element_pairs.columns = ["gene_id", "element", "value"]
    gene_element_pairs.query("value>0").to_csv("gene_element_pairs.csv", index=False)
