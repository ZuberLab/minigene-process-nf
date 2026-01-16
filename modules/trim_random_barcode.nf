process TRIM_RANDOM_BARCODE {
    tag { lane }

    input:
    tuple val(lane), path(fastq_files)

    output:
    tuple val(lane), path("output/${lane}*.fastq.gz"), emit: trimmed

    script:
    """
    barcode=\$(printf "%${params.barcode_length}s" | tr ' ' "N")
    barcode_spacer="\${barcode}${params.spacer_seq}"
    length_barcode_spacer=\${#barcode_spacer}

    mkdir -p output

    #get the read length from fastq and save it in read_length
    read_length=\$(zcat ${lane}.fastq.gz | head -n 2 | tail -n 1 | wc -c)
    read_length=\$(( read_length - 1 )) # subtract 1 for the newline character

    # calculate the minimum length of the read after trimming, subtract 4 (maximum stagger length) and the random barcode
    min_length_read=\$(( read_length - 4 - ${params.random_barcode_length} ))

    cutadapt \
        -O \${length_barcode_spacer} \
        -e ${params.spacer_error_rate} \
        -m \$min_length_read \
        -g \${barcode_spacer} \
        --action=retain \
        -j ${task.cpus} \
        -o output/${lane}.fastq.gz \
        ${lane}.fastq.gz

    """
}
