process cleanser {
        input:
            path mudata_input
            path mudata_output
            val method
            val threshold
        output:
            path mudata_output

        script:
            def thresh_opt = threshold ? "-t ${threshold}" : ""
           """
              python ${baseDir}/bin/igvf_guide_assignment.py  -i ${mudata_input} -o ${mudata_output} ${thresh_opt} --${method}
           """
}

process sceptre {
    input:
        path mudata_input
        path mudata_output
    output:
        path mudata_output

    script:
        """
            Rscript ${baseDir}/bin/sceptre_grna_assignment.R ${mudata_input} ${mudata_output}
        """
}

workflow guide_assignment {
    take:
    input_file
    output_file

    main:
    def assignment_method = params.get("ASSIGNMENT_METHOD", "cleanser").toLowerCase()

    if(assignment_method == "cleanser" || assignment_method == "umi-threshold"){
        def threshold = params.get("ASSIGNMENT_THRESHOLD", false)
        assignments = cleanser(input_file, output_file, assignment_method, threshold)
    } else if (assignment_method == "sceptre") {
        assignments = sceptre(input_file, output_file)
    }

    emit:
    assignments
}

workflow {
    guide_assignment(params.INPUT_FILE, params.OUTPUT_FILE)
}