#!/usr/bin/env nextflow

params.reference    = "/home/alexandr/Documents/octopus/data/references/genome.fna"
params.reads        = "/home/alexandr/Documents/octopus/data/reads/raw/*R{1,2}.fastq*"
params.outdir       = "results"
params.help = false

log.info """\
    ===================================
     S N P   P I P E L I N E
    ===================================
    reference: ${params.reference}
    reads    : ${params.reads}
    outdir   : ${params.outdir}
    """
    .stripIndent(true)

    
/*
 * Define the `REFINDEX` process that creates an index
 * given the reference genome file
 */
process REFINDEX {
    tag "$reference"
    input:
    path reference

    output:
    path "*"

    script:
    """
    bwa index $reference
    """
}

/*
 * Define the `QCONTROL` process that performs quality trimming and filtering of reads
 */
process QCONTROL{
    tag "${sid}"
    cpus params.cpus
    publishDir "${params.outdir}/fastp"

    input:
    tuple val(sid), path(reads)

    output:
    tuple val(sid), path(fq_1_trimmed), path(fq_2_trimmed), emit: trimmed_reads
            file("${sid}.fastp_stats.html")

    script:
    fq_1_trimmed = sid + '_R1_P.fastq.gz'
    fq_2_trimmed = sid + '_R2_P.fastq.gz'
    """
    fastp -q 20 -l 140 --trim_poly_g --thread ${task.cpus} \
    --in1 ${reads[0]} \
    --in2 ${reads[1]}\
    --out1 $fq_1_trimmed \
    --out2 $fq_2_trimmed \
    --html ${sid}.fastp_stats.html
    """
}

/*
 * Define the `ALIGN` process that aligns reads to the reference genome
 */
process ALIGN {
    tag "$reference ${sid}"
    cpus params.cpus
    publishDir "${params.outdir}/ALIGN"

    input:
    tuple val(sid), path(reads1), path(reads2)
    path reference
    path idx
    
    output:
    path "${sid}.bam"
    script:
    """
    bwa mem -t ${task.cpus} ${reference} ${reads1} ${reads2} > ${sid}.bam
    """
}

process PREPARE {
    tag "$bamFile"
	
    input:
    path bamFile

    output:
    file '*.sorted.bam'

    script:
    """
	samtools sort $bamFile -o ${bamFile.baseName}.sorted.bam
	samtools index ${bamFile.baseName}.sorted.bam
    """
}

process VARCALL {
    tag "$reference $bamFile"
    publishDir "${params.outdir}/bcftools"
	debug true
	
    input:
    path reference
    path bamFile

    output:
    file '*.vcf'

    script:
    """    
    octopus \
     -R $reference \
     -I $bamFile \
     -T chrM \
#     --config /opt/octopus/resources/configs/UMI.config \
     --sequence-error-model PCR \
     -o ${bamFile.baseName}.vcf.gz \
     --threads 8
    """
}

input_fastqs = params.reads ? Channel.fromFilePairs(params.reads, checkIfExists: true) : null

workflow {
		REFINDEX(params.reference)
		QCONTROL(input_fastqs)
       	ALIGN(QCONTROL.out[0], params.reference, REFINDEX.out)
       	PREPARE(ALIGN.out)
       	VARCALL(params.reference, PREPARE.out)
}

workflow.onComplete {
    log.info """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()
        
    log.info ( workflow.success ? "\nDone" : "\nOops" )
}




