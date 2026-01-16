process MINIGENE_MATCH {
    tag  "$sample_id"

    publishDir "${params.outputDir}/minigene_match/${library.baseName}", mode: 'copy', overwrite: true

    input:
    tuple val(lane), val(sample_id), path(fastq), path(library)

    output:
    path("${sample_id}.per_minigene_counts.tsv"), emit: minigene_counts
    path("${sample_id}.summary.tsv"), emit: minigene_counts_summary

    script:
    """
    minigene_match.py \
        --fastq ${fastq} \
        --design ${library} \
        --out-prefix ${sample_id} \
        --max-bc-mismatches 2 \
        --max-mg-mismatches ${params.minigene_mismatches} \
        --stop-len 3 \
        --start-len 3
    """
}
