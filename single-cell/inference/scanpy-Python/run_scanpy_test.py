import mudata as md
import pandas as pd
import scanpy as sc


def _run_scanpy_test_for_target(
    mdata, intended_target_name, gene_ids=None, method="wilcoxon"
):
    targeting_guides = (
        mdata["guide"].var["intended_target_name"] == intended_target_name
    )
    targeted_cells = (
        mdata["guide"][:, targeting_guides].layers["guide_assignment"].sum(axis=1)
    )
    if gene_ids is None:
        adata = mdata["gene"].copy()
    else:
        adata = mdata["gene"][:, gene_ids].copy()
    adata.obs = adata.obs.assign(
        guide_status=["present" if c > 0 else "not_present" for c in targeted_cells]
    )
    # sc.pp.log1p(mdata["gene"])
    # sc.pp.normalize_total(mdata["gene"])
    # sc.pp.regress_out(mdata["gene"], "prep_batch")
    sc.tl.rank_genes_groups(
        adata,
        groupby="guide_status",
        groups=["present"],
        method=method,
    )
    test_results_obj = adata.uns["rank_genes_groups"]
    test_results_df = pd.DataFrame(
        {
            "gene_id": test_results_obj["names"]["present"],
            "intended_target_name": intended_target_name,
            "p_value": test_results_obj["pvals"]["present"],
            "log2_fc": test_results_obj["logfoldchanges"]["present"],
        }
    )
    return test_results_df


def _run_scanpy_test_for_pairs(mdata, method="wilcoxon"):
    pairs_to_test = pd.DataFrame(mdata.uns["pairs_to_test"])
    targets = pairs_to_test["intended_target_name"].unique()
    test_results_list = []
    for target in targets:
        tested_genes = pairs_to_test.query(f"intended_target_name=='{target}'")
        test_results_target = _run_scanpy_test_for_target(
            mdata, target, tested_genes["gene_id"], method=method
        )
        test_results_list.append(test_results_target)
    mdata.uns["test_results"] = pd.concat(test_results_list).merge(
        pairs_to_test, how="right"
    )
    return mdata


def run_scanpy_test(input_mudata_fp, output_mudata_fp, method="wilcoxon"):
    mdata = md.read(input_mudata_fp)
    mdata_out = _run_scanpy_test_for_pairs(mdata, method=method)
    mdata_out.write(output_mudata_fp)
    return mdata_out
