process COMBINE_MINIGENE_MATCH_SUMMARY {
    tag { 'summarize' }

    publishDir path: "${params.outputDir}/minigene_match/",
               mode: 'copy',
               overwrite: true

    input:
    path(summaries)

    output:
    path("minigene_match_summary_all.txt"), emit: combined_summaries
    path("pair_valid_saturation.png"), emit: saturation_plot

    script:
    """
    combine_minigene_match_summary.py --plot pair_valid_saturation.png ${summaries} \
        > minigene_match_summary_all.txt
    """
}
