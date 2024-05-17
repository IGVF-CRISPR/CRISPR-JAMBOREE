nextflow.enable.dsl=2

params.mudata_input = "mudata.h5mu"
params.mudata_output = "mudata_out.h5mu"

include { guide_assignment_sceptre } from './processes/guide_assignment_sceptre.nf'

workflow {
    mudata_in = file(params.mudata_input)
    mudata_out = file(params.mudata_output)
    guide_assignment_sceptre(mudata_in, mudata_out)
}