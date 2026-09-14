nextflow.enable.dsl = 2

include { INDEX_VCF               } from '../modules/index_vcf'
include { INDEX_REFERENCE         } from '../modules/index_reference'
include { LIST_VCF_SAMPLES        } from '../modules/list_vcf_samples'
include { BCFTOOLS_EXTRACT_SAMPLE } from '../modules/bcftools_extract_sample'
include { MAP_AND_PHASE           } from '../modules/map_and_phase'
include { REBUILD_MULTISAMPLE_VCF } from '../modules/rebuild_multisample_vcf'

workflow PHASE_VCF {

    take:
    multisample_vcf
    reference
    ch_reads

    main:

    /*
     * -------------------------------------------------------------------------
     * Reference and vcf indexing
     * -------------------------------------------------------------------------
     */
    INDEX_REFERENCE(reference)
    INDEX_VCF(multisample_vcf)

    /*
    * -------------------------------------------------------------------------
    * Extract VCF sample names
    * -------------------------------------------------------------------------
    */

    LIST_VCF_SAMPLES(INDEX_VCF.out.vcf)

    /*
    * Convert sample list into meta objects
    */
    vcf_samples = LIST_VCF_SAMPLES.out.samples
        .splitText()
        .map { it.trim() }
        .filter { it }
        .map { sample ->
            [ id: sample ]
        }

    /*
    * Read samples already available from the reads channel
    */
    read_samples = ch_reads.map { meta, reads ->
        tuple(
            meta.id,
            meta,
            reads
        )
    }

    /*
    * Keep only samples present in both:
    *   - VCF
    *   - reads directory
    */
    phaseable_samples = vcf_samples
        .map { meta ->
            tuple(meta.id, meta)
        }
        .join(
            read_samples,
            by: 0
        )
        .map { sample_id, vcf_meta, read_meta, reads ->
            tuple(read_meta, reads)
        }

    /*
     * -------------------------------------------------------------------------
     * Create channel used for VCF extraction
     * -------------------------------------------------------------------------
     */
    sample_vcf_input = phaseable_samples
        .combine(INDEX_VCF.out.vcf)
        .map { meta_reads, reads, meta_vcf, vcf, tbi ->

            tuple(
                meta_reads,
                vcf,
                tbi
            )
        }

    /*
     * -------------------------------------------------------------------------
     * Extract single-sample VCFs
     * -------------------------------------------------------------------------
     */
    BCFTOOLS_EXTRACT_SAMPLE(sample_vcf_input)

    /*
     * -------------------------------------------------------------------------
     * Join extracted VCFs with ONT reads
     * -------------------------------------------------------------------------
     */
    phase_input = BCFTOOLS_EXTRACT_SAMPLE.out.vcf
        .join(phaseable_samples, by: 0)
        .map { meta, vcf, tbi, reads ->
            tuple(
                meta,
                vcf,
                tbi,
                reads
            )
        }

    /*
     * -------------------------------------------------------------------------
     * Mapping + WhatsHap phasing
     * -------------------------------------------------------------------------
     */
    reference_files = INDEX_REFERENCE.out.reference.first()
    MAP_AND_PHASE(
        reference_files,
        phase_input,
    )

    /*
     * -------------------------------------------------------------------------
     * Rebuild final multisample VCF
     * -------------------------------------------------------------------------
     */
    phased_vcf_files = MAP_AND_PHASE.out.phased_vcf
        .map { meta, vcf, tbi -> vcf }
        .collect()

    REBUILD_MULTISAMPLE_VCF(
        INDEX_VCF.out.vcf
            .map { meta, vcf, tbi ->
                tuple(vcf, tbi)
            }, 
        phased_vcfs
    )

    emit:
    phased_multisample_vcf = REBUILD_MULTISAMPLE_VCF.out.phased_vcf
}