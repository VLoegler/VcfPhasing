process INDEX_VCF {

    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta),
          path(vcf)

    output:
    tuple val(meta),
          path("${vcf.getName()}"),
          path("${vcf.getName()}.tbi"),
          emit: vcf

    script:
    """
    tabix -f -p vcf ${vcf}
    """

    stub:
    """
    touch ${vcf.getName()}.tbi
    """
}