
process starSoloGuide {
    conda "${moduleDir}/conda_envs/starsolo.yaml"

    cache 'lenient'
    cpus 2


    input:
        val WHITELIST_IN
        val READ1_GUIDE_IN
        val READ2_GUIDE_IN
        val REF_PATH_IN
    output:
        path ("Solo.out"), emit: starSolo_ouput_dir
    
    script:
        """
            STAR --soloType CB_UMI_Simple \
            --readFilesIn $READ2_GUIDE_IN $READ1_GUIDE_IN  \
            --genomeDir $REF_PATH_IN  \
            --soloCBwhitelist $WHITELIST_IN  \
            --soloCBstart 1 \
            --soloCBlen 16 \
            --soloUMIstart 17 \
            --soloUMIlen 12  \
            --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
            --outFilterScoreMin 15 \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes CR UR CY UY CB UB \
            --runThreadN ${task.cpus} \
            --soloBarcodeReadLength 0  \
            --outFilterMatchNminOverLread 0  \
            --outFilterScoreMinOverLread 0  \
            --seedSearchStartLmax 23  \
            --seedMultimapNmax 100  \
            --alignMatesGapMax 1  \
            --alignIntronMax 0  \
            --outSAMunmapped Within \
            --clipAdapterType CellRanger4
            
            
    """
}

