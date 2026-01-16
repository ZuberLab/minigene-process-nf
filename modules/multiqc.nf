process MULTIQC {
    tag { 'all' }   

    publishDir path: "${params.outputDir}",
               mode: 'copy',
               overwrite: true

    input:
    path (files)

    output:
    path("*multiqc_report.html"), emit: multiqc_report

    script:
    """
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    multiqc --interactive -f -x *.run .
    """
}