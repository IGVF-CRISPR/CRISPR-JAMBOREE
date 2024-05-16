
process soloOutGuideParsing {
        conda "${moduleDir}/conda_envs/starsolo.yaml"
        input: 
            val SOLO_OUT_DIR_IN
            val METADATA_IN 
        output:
            path ("guides_star_solo.h5ad"), emit: guides_anndata_starsolo
        script:
           """
           soloOutGuideParsing.py $SOLO_OUT_DIR_IN $METADATA_IN
           """
}

