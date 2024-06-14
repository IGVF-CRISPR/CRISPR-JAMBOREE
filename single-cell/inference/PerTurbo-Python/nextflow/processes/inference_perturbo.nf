process run_perturbo {
    input:
    path mdata

    output:
    path "mdata_out.h5mu"

    script:
    """
    run_perturbo.py -f $mdata -o mdata_out.h5mu -r $params.rna_modality -g $params.guide_modality
    """
}
