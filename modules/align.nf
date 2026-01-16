process ALIGN {
    tag { "${id}" }

    input:
    tuple val(lane), val(id), path(fastq_file), path(index)

    output:
    tuple val(lane), val(id), path("${id}.sam"), emit: alignedFiles
    path("${id}.log"), emit: alignResults
    path("${id}.sam.stats"), emit: alignStats
    path("${id}.sam.flagstat"), emit: alignFlagstats

    script:
    """
    cutadapt \
        -j ${task.cpus} \
        -l ${params.minigene_barcode_length} \
        -o minigene_barcodes_${id}.fastq.gz \
        ${fastq_file}


    bowtie2 \
        --threads \$((${task.cpus})) \
        -x ${index}/index \
        -L ${params.minigene_barcode_length} \
        --end-to-end \
        -N 0 \
        --seed 42 \
        --nofw \
        -k 1 \
        <(zcat minigene_barcodes_${id}.fastq.gz) 2> ${id}.log > ${id}.sam


    # Generate statistics for the aligned SAM file
    samtools stats ${id}.sam > ${id}.sam.stats
    samtools flagstats ${id}.sam > ${id}.sam.flagstat

    """
}
