process PROCESS_BARCODES {
    tag { barcodes.baseName }

    input:
    path(barcodes)

    output:
    path("*.fasta"), emit: processed_barcodes

    script:
    """
    process_barcodes.R ${barcodes}
    """
}
