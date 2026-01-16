process FASTQC {
    tag { id }

    input:
    tuple val(id), path(fastq_files)

    output:
    path "*_fastqc.{zip,html}", emit: fastqcResults

    script:
    """
    fastqc -t ${task.cpus} -q *.fastq.gz
    """
}