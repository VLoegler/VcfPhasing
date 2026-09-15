// modules/merge_vcfs.nf

process MERGE_VCFS {

    label "process_high"

    input:
    path vcfs
    path tbis
    path samples_order

    output:
    tuple path("multisample.phased.vcf.gz"),
          path("multisample.phased.vcf.gz.tbi"),
          emit: phased_vcf

    script:

    def vcf_list = vcfs.join(' ')

    """
    bcftools merge \
        --threads ${task.cpus} \
        ${vcf_list} | \
    bcftools view \
        --threads ${task.cpus} \
        -S ${samples_order} \
        -Oz \
        -o multisample.phased.vcf.gz

    tabix \
        -f \
        -p vcf \
        multisample.phased.vcf.gz
    """

    stub:
    """
    touch multisample.phased.vcf.gz
    touch multisample.phased.vcf.gz.tbi
    """
}