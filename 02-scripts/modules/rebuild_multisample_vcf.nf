process REBUILD_MULTISAMPLE_VCF {

    label "process_high"

    input:
    tuple path(original_vcf),
          path(original_tbi)

    path phased_vcfs

    output:
    tuple path("multisample.phased.vcf.gz"),
          path("multisample.phased.vcf.gz.tbi"),
          emit: phased_vcf

    script:

    def phased_list = phased_vcfs.join(' ')

    """
    rebuild_multisample_vcf.py \
        --input ${original_vcf} \
        --output multisample.phased.vcf.gz \
        --phased-vcfs ${phased_list}

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