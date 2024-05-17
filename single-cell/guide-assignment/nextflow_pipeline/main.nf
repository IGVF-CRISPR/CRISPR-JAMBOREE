include { guide_assignment_cleanser } from '../cleanser/nextflow/guide_assignment_cleanser.nf'
include { guide_assignment_sceptre } from '../sceptre/nextflow/processes/guide_assignment_sceptre.nf'


workflow guide_assignment {
    take:
    input_file
    output_file

    main:
    def assignment_method = params.get("ASSIGNMENT_METHOD", "cleanser").toLowerCase()

    if(assignment_method == "cleanser" || assignment_method == "umi-threshold"){
        def threshold = params.get("ASSIGNMENT_THRESHOLD", false)
        assignments = guide_assignment_cleanser(input_file, output_file, assignment_method, threshold)
    } else if (assignment_method == "sceptre") {
        assignments = guide_assignment_sceptre(input_file, output_file)
    }

    emit:
    assignments
}

workflow {
    guide_assignment(params.INPUT_FILE, params.OUTPUT_FILE)
}
