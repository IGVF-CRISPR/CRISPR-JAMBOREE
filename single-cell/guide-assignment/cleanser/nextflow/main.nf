include { guide_assignment_cleanser } from './processes/guide_assignment_cleanser.nf'


workflow guide_assignment {
    take:
    input_file
    output_file

    main:
    def assignment_method = params.get("ASSIGNMENT_METHOD", "cleanser").toLowerCase()

    def threshold = params.get("ASSIGNMENT_THRESHOLD", false)
    assignments = guide_assignment_cleanser(input_file, output_file, assignment_method, threshold)

    emit:
    assignments
}

workflow {
    guide_assignment(params.INPUT_FILE, params.OUTPUT_FILE)
}
