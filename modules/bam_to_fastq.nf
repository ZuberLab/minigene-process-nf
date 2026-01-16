process BAM_TO_FASTQ {
    tag { lane }

    input:
    tuple val(lane), path(bam)

    output:
    tuple val(lane), path("*.fastq.gz"), emit: fastq

    script:
    """
    samtools fastq -@ ${task.cpus} ${bam} -1 ${lane}_R1.fastq.gz -2 read_pairs_not_used.fastq.gz
    rm read_pairs_not_used.fastq.gz
    """
}
