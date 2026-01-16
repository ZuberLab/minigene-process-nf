process COUNT {
    tag { lane }

    publishDir path: "${params.outputDir}/counts/${saf.baseName}/${lane}",
               mode: 'copy',
               overwrite: true

    input:
    tuple val(lane), path(sams), path(saf)

    output:
    path("${lane}.txt"), emit: countedFiles
    path("${lane}.txt.summary"), emit: featureCountsResults

    script:
    """
    featureCounts \
        -T ${task.cpus} \
        -a ${saf} \
        -F SAF \
        -o ${lane}.txt \
        ${sams}
    """
}
