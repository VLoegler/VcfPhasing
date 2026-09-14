process MAP_AND_PHASE {

    tag "${meta.id}"

    label "process_high"

    input:
    tuple path(fasta),
          path(fai),
          path(mmi)

    tuple val(meta),
          path(vcf),
          path(tbi),
          path(reads)

    output:
    tuple val(meta),
          path("${meta.id}.phased.vcf.gz"),
          path("${meta.id}.phased.vcf.gz.tbi"),
          emit: phased_vcf

    script:

    def sort_mem = (task.memory.giga / task.cpus / 2).intValue()

    """
    minimap2 \
        -t ${task.cpus} \
        -ax map-ont \
        ${mmi} \
        ${reads} \
    | samtools sort \
        -@ ${task.cpus} \
        -m ${sort_mem}G \
        -o ${meta.id}.bam

    samtools index \
        -@ ${task.cpus} \
        ${meta.id}.bam

    whatshap phase \
        --reference ${fasta} \
        --output ${meta.id}.phased.vcf.gz \
        --ignore-read-groups \
        --indels \
        ${vcf} \
        ${meta.id}.bam

    tabix -f -p vcf ${meta.id}.phased.vcf.gz

    rm -f ${meta.id}.bam
    rm -f ${meta.id}.bam.bai
    """

    stub:
    """
    touch ${meta.id}.phased.vcf.gz
    touch ${meta.id}.phased.vcf.gz.tbi
    """
}