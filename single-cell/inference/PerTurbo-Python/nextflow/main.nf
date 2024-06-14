#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.mdata_fp = "$baseDir/data/gasperini_pilot_subset.h5mu"

include { run_perturbo } from './processes/inference_perturbo.nf'

workflow {
    Channel.fromPath(params.mdata_fp) | run_perturbo
}