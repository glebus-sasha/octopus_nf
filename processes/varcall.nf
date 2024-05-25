// Define the `VARCALL` process that performs variant calling
process VARCALL {
    container = 'dancooke/octopus:sandybridge'
    tag "$reference $bamFile"
    publishDir "${params.outdir}/octopus"
	debug true
//    errorStrategy 'ignore'
	
    input:
    path reference
    tuple val(sid), path(bai), path(bamFile)
    path fai

    output:
    tuple val(sid), path("${sid}.vcf.gz"),      emit:vcf

    script:
    """    
    octopus \
    -R $reference \
    -I $bamFile \
    --sequence-error-model PCR \
    --forest /opt/octopus/resources/forests/germline.v0.7.4.forest \
    -o ${sid}.vcf.gz \
    --threads ${task.cpus}  \
    -C cancer \
    --very-fast
    """
}