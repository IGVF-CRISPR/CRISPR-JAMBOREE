
process parsing_guide_metadata {
        conda "${moduleDir}/conda_envs/starsolo.yaml"
        input: 
            val guide_metadata
        output:
            path ("guides.bed"), emit: guides_bed
            path ("isoforms.txt"), emit: pseudo_isoforms_guides
            path ("guides.gtf"), emit: guides_gtf
            path ("pseudo_genome.fa"), emit : pseudo_genome
        script:
           """
              parsing_guide_medata.py  $guide_metadata
              bed2gtf --bed guides.bed --output guides.gtf --isoforms 'isoforms.txt'
           """
}
