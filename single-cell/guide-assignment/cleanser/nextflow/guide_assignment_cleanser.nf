process guide_assignment_cleanser {
    input:
        path mudata_input
        path mudata_output
        val method
        val threshold
    output:
        path mudata_output, emit: mudata_output

    script:
        def thresh_opt = threshold ? "-t ${threshold}" : ""
        """
            python ${moduleDir}/bin/igvf_guide_assignment.py  -i ${mudata_input} -o ${mudata_output} ${thresh_opt} --${method}
        """
}
