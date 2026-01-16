#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def helpMessage() {
    log.info"""
    ================================================================
     dual-crispr-process-nf
    ================================================================
     DESCRIPTION

     Process CRISPR and shRNA functional genetic screening data.

     Usage:
     nextflow run zuberlab/crispr-process-nf

     Options:
        --inputDir                          Input directory containing raw files. Either BAM or FASTQ files.
                                            The FASTQ files must be named <lane>.fastq.gz.
                                            The BAM files must be named <lane>.bam.
                                            (default: '01_raw')

        --outputDir                         Output directory for processed files.
                                            (default: '02_processed')

        --pheno                             Path to metadata file (default pheno.xslx)

        --library                           Path to minigene library file.
                                            (default: 'library.txt')
                                            The following columns are required:
                                                - id:       unique name of sgRNA
                                                - gene:     gene targeted by sgRNA
                                                - sequence: nucleotide sequence of sgRNA

         --barcodes                         Path to file containing barcodes for demultiplexing.
                                            (default: 'barcodes.fasta')
                                            The following columns are required:
                                                - lane:         name of BAM / FASTQ input file
                                                - sample_name:  name of demultiplexed sample
                                                - barcode:   nucleotide sequence of the sample barcode


        --barcode_demux_mismatches          Number of mismatches allowed during demultiplexing
                                            of barcode. (default: 1)


        --barcode_length                    Number of nucleotides in sample barcode.
                                            (default: 4)

        --spacer_seq                        Nucleotide sequence of spacer between
                                            barcodes and minigene sequence.
                                            (default: AGTATTAGGCGTCAAGGTCCTTA)

        --spacer_error_rate                 Error rate for spacer sequence. (default: 3)

        --minigene_barcode_length           Number of nucleotides in minigene barcode. (default: 8)

        --minigene_length                   Number of nucleotides in minigene sequence. (default: 99)

        --minigene_mismatches               Number of allowed mismatches in minigene sequence. (default: 9)

        --random_barcode_length             Length of random barcode. (default: 4)

     Profiles:
        standard                    local execution
        apptainer                   local execution with apptainer
        cbe                         SLURM execution with apptainer on CBE cluster


     Docker:
     zuberlab/dual-crispr-nf:1.1

     Author:
     Florian Andersch (florian.andersch@imp.ac.at)


    """.stripIndent()
}

params.help = false

if (params.help) {
    helpMessage()
    exit 0
}

log.info ""
log.info " parameters "
log.info " =============================================================="
log.info " input directory                                  : ${params.inputDir}"
log.info " output directory                                 : ${params.outputDir}"
log.info " metadata file                                    : ${params.pheno}"
log.info " library file                                     : ${params.library}"
log.info " barcode file                                     : ${params.barcodes}"
log.info " barcode length                                   : ${params.barcode_length}"
log.info " mismatch allowance for demultiplexing            : ${params.barcode_demux_mismatches}"
log.info " spacer seq (nt)                                  : ${params.spacer_seq}"
log.info " spacer error rate                                : ${params.spacer_error_rate}"
log.info " minigene barcode length                          : ${params.minigene_barcode_length}"
log.info " minigene length (nt)                             : ${params.minigene_length}"
log.info " mismatch allowance for minigene sequence         : ${params.minigene_mismatches}"
log.info " random barcode length                            : ${params.random_barcode_length}"
log.info " =============================================================="
log.info ""

// Import modules
include { BAM_TO_FASTQ } from './modules/bam_to_fastq'
include { TRIM_RANDOM_BARCODE } from './modules/trim_random_barcode'
include { PROCESS_BARCODES } from './modules/process_barcodes'
include { DEMULTIPLEX } from './modules/demultiplex'
include { TRIM_BARCODE_AND_SPACER } from './modules/trim_barcode_and_spacer'
include { MINIGENE_MATCH } from './modules/minigene_match'
include { COMBINE_MINIGENE_MATCH_SUMMARY } from './modules/combine_minigene_match_summary'
include { PROCESS_LIBRARY } from './modules/process_library'
include { BOWTIE_INDEX } from './modules/bowtie_index'
include { ALIGN } from './modules/align'
include { COUNT } from './modules/count'
include { COMBINE_COUNTS } from './modules/combine_counts'
include { FASTQC } from './modules/fastqc'
include { MULTIQC } from './modules/multiqc'
include { PCA } from './modules/pca'

// Define input channels
ch_input_bam = Channel.fromPath("${params.inputDir}/*.bam")
    .map { file -> tuple(file.baseName, file) }

ch_input_fastq = Channel.fromPath("${params.inputDir}/*.fastq.gz")
    .map {
        file ->
            def lane = file.name.toString().replaceAll(/\.fastq\.gz$/, '')
            tuple(lane, file)
        }
    .groupTuple()

ch_barcodes = Channel.fromPath(params.barcodes)
ch_library = Channel.fromPath(params.library)
ch_pheno = Channel.fromPath(params.pheno)

// Main workflow
workflow {
    // BAM to FASTQ conversion
    ch_fastq_from_bam = BAM_TO_FASTQ(ch_input_bam)

    // Combine BAM-derived and direct FASTQ inputs
    ch_all_fastq = ch_fastq_from_bam.mix(ch_input_fastq)

    // Trim random barcode
    ch_trimmed_random = TRIM_RANDOM_BARCODE(ch_all_fastq)

    // Process barcodes
    ch_processed_barcodes = PROCESS_BARCODES(ch_barcodes)

    // Combine processed barcodes with trimmed random barcodes
    ch_trimmed_random_barcodes = ch_processed_barcodes
        .flatten()
        .map { barcode ->
            def lane = barcode.name.toString().replaceAll(/\.fasta$/, '')
            [lane, barcode]
        }
        .groupTuple()
        .join(ch_trimmed_random)

    // Demultiplex
    ch_demuxed = DEMULTIPLEX(ch_trimmed_random_barcodes)

    // Flatten demultiplexed files
    ch_demuxed_flattened = ch_demuxed
        .flatMap { lane, files ->
                files.collect { file ->
                    def id = file.name.toString().replaceAll(/\.fastq\.gz$/, '')
                    [lane, id, file]
                }
        }.groupTuple(by: [0,1])

    // Split the channel based on whether id ends with "unknown"
    ch_demuxed_flattened
        .branch {
            unknown: it[1].toString().endsWith("#unknown")
            known: true
        }
        .set { ch_demuxed_flattened_split }

    // Trim barcode and spacer
    ch_trimmed_spacer = TRIM_BARCODE_AND_SPACER(ch_demuxed_flattened_split.known)

    ch_trimmed_spacer_combined_minigene = ch_trimmed_spacer
        .combine(ch_library)

    // match minigenes
    ch_minigene_matched = MINIGENE_MATCH(ch_trimmed_spacer_combined_minigene)

    // Combine minigene match stats
    COMBINE_MINIGENE_MATCH_SUMMARY(ch_minigene_matched.minigene_counts_summary.collect())

    // Process library
    ch_library_out = PROCESS_LIBRARY(ch_library)

    // Bowtie index
    ch_bt2_index = BOWTIE_INDEX(ch_library_out.fasta)

    // Combine trimmed spacer files and bowtie index
    ch_trimmed_spacer_combined_index = ch_trimmed_spacer
        .combine(ch_bt2_index)

    // Align
    ch_aligned = ALIGN(ch_trimmed_spacer_combined_index)

    // Group aligned files by lane
    ch_grouped_aligned =  ch_aligned.alignedFiles
        .map { lane, id, file -> tuple(lane, file) }
        .groupTuple()
        .combine(ch_library_out.saf)

    // Count
    ch_counted = COUNT(ch_grouped_aligned)

    // Combine counts
    ch_counts_pca = COMBINE_COUNTS(ch_counted.countedFiles.collect(), ch_library)

    //PCA
    if (params.pheno && file(params.pheno).exists()) {
        PCA(ch_counts_pca.combine(ch_pheno))
    } else {
        log.info "Skipping PCA: No phenotype file provided"
    }

    // collect all fastq files
    ch_fastq_files = ch_all_fastq
        .mix(ch_demuxed_flattened.map { lane, baseName, file -> tuple(baseName, file) })

    // FastQC
    ch_fastqc = FASTQC(ch_fastq_files)

    // Combine FASTQC files, alignment results, and featureCounts results
    ch_multiqqc_files = ch_fastqc.collect()
        .mix(ch_aligned.alignResults.collect())
        .mix(ch_aligned.alignStats.collect())
        .mix(ch_aligned.alignFlagstats.collect())
        .mix(ch_counted.featureCountsResults.collect())
        .collect()

    // MultiQC
    MULTIQC(ch_multiqqc_files)
}

// On completion
workflow.onComplete {
    println ( workflow.success ? "COMPLETED!" : "FAILED" )
}
