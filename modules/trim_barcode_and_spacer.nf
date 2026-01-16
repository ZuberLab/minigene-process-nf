process TRIM_BARCODE_AND_SPACER {
    tag { id }

    input:
    tuple val(lane), val(id), path(fastq_file)

    output:
    tuple val(lane), val(id), path("output/${lane}*.fastq.gz"), emit: trimmed

    script:
    """
    barcode=\$(printf "%${params.barcode_length}s" | tr ' ' "N")
    barcode_spacer="\${barcode}${params.spacer_seq}"
    length_barcode_spacer=\${#barcode_spacer}
    length_minigene=\$(( ${params.minigene_length} + ${params.minigene_barcode_length} ))

    mkdir -p output

    cutadapt \
        -j ${task.cpus} \
        -u \${length_barcode_spacer} \
        -l \${length_minigene} \
        --minimum-length \${length_minigene} \
        -o output/${id}.fastq.gz \
        ${fastq_file}
    """
}
