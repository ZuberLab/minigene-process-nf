process PCA {
    tag { 'all' }

    publishDir path: "${params.outputDir}/",
               mode: 'copy',
               overwrite: true

    input:
    tuple path(counts), path(pheno)

    output:
    path("pca/*"), emit : pca_files

    script:
    """
    pca.R ${counts} ${pheno}
    """
}
