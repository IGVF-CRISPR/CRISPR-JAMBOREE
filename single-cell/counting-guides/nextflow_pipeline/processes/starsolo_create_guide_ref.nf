
process starsolo_create_guide_ref {
        conda "${moduleDir}/conda_envs/starsolo.yaml"
        input: 
            val pseudo_genome_in
            val guides_gtf_in 
        output:
            path ("new_ref"), emit: reference_guide_dir
        script:
           """
           mkdir new_ref
           STAR --runThreadN 4 --runMode genomeGenerate --genomeDir new_ref  --genomeFastaFiles $pseudo_genome_in --genomeSAindexNbases 5 --sjdbGTFfile $guides_gtf_in 
           """
}

