nextflow.enable.dsl=2

params.mudata_fp = "mudata.h5mu"

include { inference_sceptre } from './processes/inference_sceptre.nf'

workflow {
    mudata_in = file(params.mudata_fp)
    inference_sceptre(mudata_in)
}