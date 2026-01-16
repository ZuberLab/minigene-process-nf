process PROCESS_LIBRARY {
    tag { library.baseName }

    input:
    path(library)

    output:
    path("${library.baseName}.saf"), emit : saf
    path("${library.baseName}.fasta"), emit : fasta

    script:
    """
    process_library.R ${library}
    """
}
