process BCFTOOLS_EXTRACT_SAMPLE {

    tag "${meta.id}"
    label "process_low"

    input:
    tuple val(meta),
          path(vcf),
          path(tbi)

    output:
    tuple val(meta),
          path("${meta.id}.vcf.gz"),
          path("${meta.id}.vcf.gz.tbi"),
          emit: vcf

    script:
    """
    bcftools view \
        --threads ${task.cpus} \
        -s ${meta.id} \
        -Oz \
        -o ${meta.id}.vcf.gz \
        ${vcf}

    tabix -f -p vcf ${meta.id}.vcf.gz
    """

    stub:
    """
    touch ${meta.id}.vcf.gz
    touch ${meta.id}.vcf.gz.tbi
    """
}