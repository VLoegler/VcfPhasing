nextflow.enable.dsl = 2

include { INDEX_VCF               } from '../modules/index_vcf'
include { INDEX_REFERENCE         } from '../modules/index_reference'
include { LIST_VCF_SAMPLES        } from '../modules/list_vcf_samples'
include { BCFTOOLS_EXTRACT_SAMPLE } from '../modules/bcftools_extract_sample'
include { MAP_AND_PHASE           } from '../modules/map_and_phase'
include { MERGE_VCFS              } from '../modules/merge_vcfs'

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
    vcf_ids = vcf_samples
        .map { it.id }

    /*
    * Read samples already available from the reads channel
    */
    read_sample_ids = ch_reads
        .map { meta, reads ->
            tuple(
                meta.id,
                reads
            )
        }
    valid_reads = read_sample_ids
        .join(
            vcf_ids.map { id -> tuple(id, true) },
            by: 0
        )
        .map { sample_id, reads, dummy ->
            tuple(sample_id, reads)
        }

    /*
     * -------------------------------------------------------------------------
     * Create channel used for VCF extraction
     * -------------------------------------------------------------------------
     */
    sample_vcf_input = vcf_samples
        .combine(INDEX_VCF.out.vcf)
        .map { meta_sample, meta_vcf, vcf, tbi ->

            tuple(
                meta_sample,
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
    sample_vcfs = BCFTOOLS_EXTRACT_SAMPLE.out.vcf
        .map { meta, vcf, tbi ->
            tuple(
                meta.id,
                meta,
                vcf,
                tbi
            )
        }

    /*
     * -------------------------------------------------------------------------
     * Join extracted VCFs with ONT reads
     * -------------------------------------------------------------------------
     */
    phase_input = sample_vcfs
        .join(valid_reads, by:0)
        .map { sample_id, meta, vcf, tbi, reads ->

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
     * Identify unphaseable samples
     * -------------------------------------------------------------------------
     */
    sample_vcfs
        .join(
            read_sample_ids,
            by: 0,
            remainder: true
        )
        .view()
    unphaseable_vcfs = sample_vcfs
        .join(
            valid_reads,
            by:0,
            remainder:true
        )
        .filter { sample_id, meta, vcf, tbi, reads ->
            reads == null
        }
        .map { sample_id, meta, vcf, tbi, reads ->

            tuple(
                meta,
                vcf,
                tbi
            )
        }

    /*
     * -------------------------------------------------------------------------
     * Rebuild final multisample VCF
     * -------------------------------------------------------------------------
     */
    merged_vcfs =
        MAP_AND_PHASE.out.phased_vcf
        .mix(unphaseable_vcfs)
        .map { meta, vcf, tbi -> vcf }
        .collect()
    sample_order = LIST_VCF_SAMPLES.out.samples

    MERGE_VCFS(
        merged_vcfs, 
        sample_order,
    )

    emit:
    phased_multisample_vcf = MERGE_VCFS.out.phased_vcf
}