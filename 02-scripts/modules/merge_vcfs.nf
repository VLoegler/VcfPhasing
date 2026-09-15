// modules/merge_vcfs.nf

process MERGE_VCFS {

    label "process_high"

    input:
    val chr
    path vcfs
    path tbis
    path samples_order

    output:
    tuple path("${chr}.multisample.phased.vcf.gz"),
          path("${chr}.multisample.phased.vcf.gz.tbi"),
          emit: phased_vcf

    script:

    def vcf_list = vcfs.join(' ')

    """
    ulimit -n 4096
    
    bcftools merge \
        -r ${chr} \
        --threads ${task.cpus} \
        ${vcf_list} | \
    bcftools view \
        --threads ${task.cpus} \
        -S ${samples_order} \
        -Oz \
        -o ${chr}.multisample.phased.vcf.gz

    tabix \
        -f \
        -p vcf \
        ${chr}.multisample.phased.vcf.gz
    """

    stub:
    """
    touch ${chr}.multisample.phased.vcf.gz
    touch ${chr}.multisample.phased.vcf.gz.tbi
    """
}