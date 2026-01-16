process COMBINE_COUNTS {
    tag { library.baseName }

    publishDir path: "${params.outputDir}/counts/${library.baseName}",
               mode: 'copy',
               overwrite: true

    input:
    path(counts)
    path(library)

    output:
    path("counts_mageck.txt"), emit: combined_counts

    script:
    """
    combine_counts.R ${library} ${counts} > counts_mageck.txt
    """
}