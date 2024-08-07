#!/usr/bin/env nextflow
// Include processes
include { ALIGN }               from './processes/align.nf'
include { FLAGSTAT }            from './processes/flagstat.nf'
include { QUALIMAP }            from './processes/qualimap.nf'
include { BAMINDEX }            from './processes/bamindex.nf'
include { VARCALL }             from './processes/varcall.nf'
include { ANNOTATE }            from './processes/annotate.nf'
include { REPORT }              from './processes/report.nf'

// Logging pipeline information
log.info """\
\033[0;36m  ==========================================  \033[0m
\033[0;34m       O C T O P U S   P I P E L I N E        \033[0m
\033[0;36m  ==========================================  \033[0m

    reference:  ${params.reference}
    reads:      ${params.reads}
    outdir:     ${params.outdir}
    """
    .stripIndent(true)

// Define help
if ( params.help ) {
    help = """main.nf: This repository contains a Nextflow pipeline for analyzing 
            |Next-Generation Sequencing (NGS) data using octopus 
            |
            |Required arguments:
            |   --reference     Location of the reference file.
            |                   [default: ${params.reference}]
            |   --reads         Location of the input file file.
            |                   [default: ${params.reads}]
            |   --outdir        Location of the output file file.
            |                   [default: ${params.outdir}]
            |
            |Optional arguments:
            |   -profile        <docker/singularity>
            |   -prebuild       Use pre-built bwa indexes and pre-downloaded vep cache
            |   -reports        Generate pipeline reports
            |
""".stripMargin()
    // Print the help with the stripped margin and exit
    println(help)
    exit(0)
}

// Define the input channel for FASTQ files, if provided
input_fastqs = params.reads ? Channel.fromFilePairs("${params.reads}/*[rR]{1,2}*.*{fastq,fq}*", checkIfExists: true) : null

// Define the input channel for bwa index files, if provided
bwaidx = params.bwaidx ? Channel.fromPath("${params.bwaidx}/*", checkIfExists: true).collect() : null

// Define the input channel for fai index files, if provided
faidx = params.bwaidx ? Channel.fromPath("${params.faidx}/*.fai", checkIfExists: true).collect() : null


// Define the input channels for Clinvar files and indeces, if provided
//clinvar_gz = params.bwaidx ? Channel.fromPath("${params.vepcache}/clinvar.vcf.gz", checkIfExists: true) : null
//clinvar_gz_tbi = params.bwaidx ? Channel.fromPath("${params.vepcache}/clinvar.vcf.gz.tbi", checkIfExists: true) : null

// Define the bed_file channel
// If params.regions is provided, create a channel from the specified path and collect it into a list
// Otherwise, create a channel from the path "assets/dummy.bed" and collect it into a list
// https://github.com/nextflow-io/nextflow/issues/1694
bed_file = params.regions ? Channel.fromPath("${params.regions}").collect() : Channel.fromPath("assets/dummy.bed").collect()

// Define the workflow
workflow {
    ALIGN(input_fastqs, params.reference, bwaidx, bed_file)
    FLAGSTAT(ALIGN.out.bam)
    QUALIMAP(ALIGN.out.bam)
    BAMINDEX(ALIGN.out.bam)
    VARCALL(params.reference, BAMINDEX.out.bai, faidx, bed_file)
    ANNOTATE(VARCALL.out.vcf)
    REPORT(FLAGSTAT.out.flagstat.collect(), QUALIMAP.out.collect(), ANNOTATE.out.html.collect())
    
    // Make the pipeline reports directory if it needs
    if ( params.reports ) {
        def pipeline_report_dir = new File("${params.outdir}/pipeline_info/")
        pipeline_report_dir.mkdirs()
    }
}

// Log pipeline execution summary on completion
workflow.onComplete {
    log.info """\033[0;32m\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete.format('yyyy-MM-dd_HH-mm-ss')}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        \033[0m"""
        .stripIndent()
        
    log.info ( workflow.success ? "\nDone" : "\nOops" )
}