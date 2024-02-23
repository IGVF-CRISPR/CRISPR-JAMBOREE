import mudata as md
import pandas as pd
import numpy as np
import perturbo


def run_perturbo(mdata_input_fp, mdata_output_fp):
    """
    Run PerTurbo on the selected guide--element pairs and return a new MuData object
    with the test results stored in mdata.uns["test_results"]
    """
    mdata = md.read(mdata_input_fp)
    mdata["gene"].obs = (
        mdata.obs.join(mdata["gene"].obs)
        .join(mdata["guide"].obs)
        .assign(log1p_total_guide_umis=lambda x: np.log1p(x["total_guide_umis"]))
    )
    mdata["guide"].X = mdata["guide"].layers["guide_assignment"]
    pairs_to_test_df = pd.DataFrame(mdata.uns["pairs_to_test"])
    mdata.uns["intended_target_names"] = sorted(
        pd.unique(pairs_to_test_df["intended_target_name"])
    )

    mdata["gene"].varm["intended_targets"] = (
        pairs_to_test_df.assign(value=1)
        .pivot(index="gene_id", columns="intended_target_name", values="value")
        .reindex(mdata["gene"].var_names)
        .fillna(0)
    )

    mdata.uns["intended_target_names"] = sorted(
        pd.unique(pairs_to_test_df["intended_target_name"])
    )

    intended_targets_df = pd.get_dummies(
        mdata["guide"].var["intended_target_name"]
    ).astype(float)

    mdata["guide"].varm["intended_targets"] = intended_targets_df[
        mdata.uns["intended_target_names"]
    ]

    perturbo.PERTURBO.setup_mudata(
        mdata,
        batch_key="prep_batch",
        library_size_key="total_gene_umis",
        continuous_covariates_keys=["total_guide_umis"],
        guide_by_element_key="intended_targets",
        gene_by_element_key="intended_targets",
        modalities={
            "rna_layer": "gene",
            "perturbation_layer": "guide",
        },
    )

    model = perturbo.PERTURBO(mdata, likelihood="nb")
    model.train(20, lr=0.01, batch_size=128)

    igvf_name_map = {
        "element": "intended_target_name",
        "gene": "gene_id",
        "q_value": "posterior_probability",
    }

    element_effects = (
        model.get_element_effects()
        .rename(columns=igvf_name_map)
        .assign(log2_fc=lambda x: x["loc"] / np.log(2))
        .merge(pairs_to_test_df)
    )

    mdata = md.read(mdata_input_fp)
    mdata.uns["test_results"] = element_effects[
        [
            "intended_target_name",
            "gene_id",
            "posterior_probability",
            "log2_fc",
        ]
    ]
    mdata.write(mdata_output_fp)
    return mdata
