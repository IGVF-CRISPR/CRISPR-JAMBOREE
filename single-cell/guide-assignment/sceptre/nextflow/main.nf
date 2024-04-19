nextflow.enable.dsl=2

params.mudata_fp = "mudata.h5mu"

include { guide_assignment_sceptre } from './processes/guide_assignment_sceptre.nf'

workflow {
    mudata_in = file(params.mudata_fp)
    guide_assignment_sceptre(mudata_in)
}