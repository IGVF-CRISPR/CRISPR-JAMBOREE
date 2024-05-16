process inference_sceptre {
    container 'igvf/sceptre-igvf:v0.2'

    input:
    path mudata_fp

    output:
    path "mudata_out.h5mu"

    script:
    """
    inference_sceptre.R \
      $mudata_fp \
      $params.side \
      $params.grna_integration_strategy \
      $params.resampling_approximation \
      $params.control_group \
      $params.resampling_mechanism \
      '${params.formula_object}'
    """
}
