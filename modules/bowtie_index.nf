process BOWTIE_INDEX {
    tag { library_fasta.baseName }

    input:
    path(library_fasta)

    output:
    path("bt2"), emit: bt2Index

    script:
    """
    mkdir -p bt2
    bowtie2-build ${library_fasta} bt2/index
    """
}