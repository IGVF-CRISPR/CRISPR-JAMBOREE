nextflow.enable.dsl=2

include { parsing_guide_metadata }    from  './processes/parsing_guide_metadata'
include { starsolo_create_guide_ref } from './processes/starsolo_create_guide_ref'
include { starSoloGuide }             from './processes/starSoloGuide'
include { soloOutGuideParsing }       from './processes/soloOutGuideParsing'





workflow {
    metadata_out = parsing_guide_metadata(params.GUIDES_METADATA)
    ref_dir = starsolo_create_guide_ref(metadata_out.pseudo_genome, metadata_out.guides_gtf )
    star_guide_out = starSoloGuide(params.WHITELIST, params.FASTQ_FILES_GUIDES[0], params.FASTQ_FILES_GUIDES[1], ref_dir.reference_guide_dir)
    ann_out = soloOutGuideParsing(star_guide_out.starSolo_ouput_dir, params.GUIDES_METADATA  )

}

