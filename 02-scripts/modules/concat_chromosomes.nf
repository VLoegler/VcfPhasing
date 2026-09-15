process CONCAT_CHROMOSOMES {

    label "process_medium"

    input:
    path vcfs
    path tbis

    output:
    tuple path("multisample.phased.vcf.gz"),
          path("multisample.phased.vcf.gz.tbi"),
          emit: vcf

    script:

    def orderedVcfs = vcfs
        .sort { a, b ->

            def chrA = (
                a.baseName =~ /chromosome(\d+)/
            )[0][1] as Integer

            def chrB = (
                b.baseName =~ /chromosome(\d+)/
            )[0][1] as Integer

            chrA <=> chrB
        }
        .join(' ')

    """
    bcftools concat \
        --threads ${task.cpus} \
        ${orderedVcfs} \
        -Oz \
        -o multisample.phased.vcf.gz

    tabix -f -p vcf multisample.phased.vcf.gz
    """
}