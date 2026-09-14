/*
 * Author: Victor Loegler
 * Date: 2026-09-11
 * Description:
 * Phase variants in a multisample VCF using ONT reads and WhatsHap.
 */

nextflow.enable.dsl = 2

include { PHASE_VCF } from './workflows/phase_vcf'

workflow {

    // Input VCF
    ch_multisample_vcf = Channel
        .fromPath(
            params.multisample_vcf,
            checkIfExists: true
        )
        .map { vcf ->
            tuple(
                [id: 'multisample'],
                vcf
            )
        }

    // Reference genome
    ref_fasta = Channel.fromPath(
        params.ref_fasta,
        checkIfExists: true
    )

    // Discover ONT read files
    ch_reads = Channel
        .fromPath(
            "${params.reads_dir}/*.fastq.gz",
            checkIfExists: true
        )
        .map { reads ->

            def sample_id =
                reads.baseName
                     .replaceFirst(/\.fastq$/, '')

            def meta = [
                id : sample_id
            ]

            tuple(meta, reads)
        }

    // Run phasing workflow
    PHASE_VCF(
        ch_multisample_vcf,
        ref_fasta,
        ch_reads,
    )
}