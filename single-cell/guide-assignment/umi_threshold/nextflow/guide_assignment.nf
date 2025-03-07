process threshold {
    container 'igvf/cleanser:v1.2'

    input:
        path mudata_input
        path mudata_output
        val threshold
    output:
        path mudata_output, emit: mudata_output

    script:
        def thresh_opt = threshold ? "-t ${threshold}" : ""
        """
            python ${moduleDir}/bin/threshold_assignment.py  -i ${mudata_input} -o ${mudata_output} ${thresh_opt}
        """
}
