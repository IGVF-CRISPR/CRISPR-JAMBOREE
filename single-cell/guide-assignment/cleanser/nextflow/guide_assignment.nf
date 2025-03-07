process cleanser {
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
            cleanser -i ${mudata_input} --posteriors-output ${mudata_output} --modality guide --capture-method capture_method --output-layer guide_assignment ${thresh_opt}
        """
}
