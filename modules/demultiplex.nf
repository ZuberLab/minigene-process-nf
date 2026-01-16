process DEMULTIPLEX {
    tag { lane }

    input:
    tuple val(lane), path(barcodes), path(fastq_file)

    output:
    tuple val(lane), path("*.fastq.gz"), emit: demuxed

    script:
    """
    cutadapt \
        -j ${task.cpus} \
        -e ${params.barcode_demux_mismatches} \
        --no-indels \
        -g file:${barcodes} \
        --action=none \
        -o "${lane}#{name}.fastq.gz" \
        ${fastq_file}
    """
}
