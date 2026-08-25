process ALIGN {
    tag { "${id}" }

    input:
    tuple val(lane), val(id), path(fastq_file), path(index)

    output:
    tuple val(lane), val(id), path("${id}.marked.sam"), emit: alignedFiles
    path("${id}.log"), emit: alignResults
    path("${id}.marked.sam.stats"), emit: alignStats
    path("${id}.marked.sam.flagstat"), emit: alignFlagstats

    script:
    """
    cutadapt -j ${task.cpus} -l ${params.minigene_alignment_length} ${fastq_file} | \
    seqkit seq -r -p | \
    bowtie2 --threads ${task.cpus} \
        -x ${index}/index \
        -L ${params.minigene_barcode_length} \
        --end-to-end \
        --ignore-quals \
        --rdg 1000,1000 \
        --rfg 1000,1000 \
        --score-min C,${params.minigene_alignment_score},0 \
        -N 0 \
        -i S,1,0 \
        --seed 42 \
        --norc \
        - 2> ${id}.log > ${id}.sam

    mark_sam_by_position_mismatches.py \
        --input ${id}.sam \
        --output ${id}.marked.sam \
        --prefix-length ${params.minigene_barcode_length} \
        --max-prefix-mismatches 1 \
        --max-suffix-mismatches 5 \
        --verbose

    # Generate statistics for the aligned SAM file
    samtools stats ${id}.marked.sam > ${id}.marked.sam.stats
    samtools flagstat ${id}.marked.sam > ${id}.marked.sam.flagstat

    """
}
