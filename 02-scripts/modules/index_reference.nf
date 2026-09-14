process INDEX_REFERENCE {

    tag "${fasta.baseName}"
    label 'process_low'

    input:
    path fasta

    output:
    tuple path("reference.fasta"),
          path("reference.fasta.fai"),
          path("reference.mmi"),
          emit: reference

    script:
    """
    cp ${fasta} reference.fasta

    samtools faidx reference.fasta

    minimap2 \
        -d reference.mmi \
        reference.fasta
    """

    stub:
    """
    touch reference.fasta
    touch reference.fasta.fai
    touch reference.mmi
    """
}