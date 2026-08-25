process PCA {
    tag { 'all' }

    publishDir path: "${params.outputDir}/pcas/${library.baseName}",
               mode: 'copy',
               overwrite: true

    input:
    tuple path(counts), path(pheno)
    path(library)


    output:
    path("pca/*"), emit : pca_files

    script:
    """
    pca.R ${counts} ${pheno}
    """
}
