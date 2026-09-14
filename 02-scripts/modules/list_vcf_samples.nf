process LIST_VCF_SAMPLES {

    label 'process_low'

    input:
    tuple val(meta),
          path(vcf),
          path(tbi)

    output:
    path("samples.txt"),
    emit: samples

    script:
    """
    bcftools query -l ${vcf} > samples.txt
    """

    stub:
    """
    touch samples.txt
    """
}