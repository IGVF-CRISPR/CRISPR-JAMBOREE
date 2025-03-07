include { cleanser ; threshold } from './guide_assignment.nf'


workflow guide_assignment {
    take:
    input_file
    output_file

    main:
    def assignment_method = params.get("ASSIGNMENT_METHOD", "cleanser").toLowerCase()
    def threshold_value = params.get("ASSIGNMENT_THRESHOLD", false)

    if( assignment_method == "cleanser" ) {
        assignments = cleanser(input_file, output_file, threshold_value)
    } else {
        assignments = threshold(input_file, output_file, threshold_value)
    }

    emit:
    assignments
}

workflow {
    guide_assignment(params.INPUT_FILE, params.OUTPUT_FILE)
}
