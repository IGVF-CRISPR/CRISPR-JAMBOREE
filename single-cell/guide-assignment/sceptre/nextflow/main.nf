nextflow.enable.dsl=2

params.mudata_fp = "mudata.h5mu"

process guide_assignment_sceptre {
    container 'ekatsevi/sceptre-igvf:v0.1'

    input:
    path mudata_fp

    output:
    path "mudata_out.h5mu"

    script:
    """
    Rscript -e "
      mudata_in <- MuData::readH5MU('${mudata_fp}')
      mudata_out <- sceptreIGVF::assign_grnas_sceptre(mudata = mudata_in)
      MuData::writeH5MU(mudata_out, 'mudata_out.h5mu')
    "
    """
}

workflow {
    mudata_in = file(params.mudata_fp)
    guide_assignment_sceptre(mudata_in)
}