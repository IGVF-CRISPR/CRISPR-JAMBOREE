process guide_assignment_sceptre {
    container 'igvf/sceptre-igvf:v0.2'

    input:
        path mudata_input
        path mudata_output
    output:
        path mudata_output

    script:
    """
    Rscript -e "
      mudata_in <- MuData::readH5MU('${mudata_input}')
      mudata_out <- sceptreIGVF::assign_grnas_sceptre(mudata = mudata_in)
      MuData::writeH5MU(mudata_out, '${mudata_output}')
    "
    """
}
