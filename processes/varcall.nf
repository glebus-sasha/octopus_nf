// Define the `VARCALL` process that performs variant calling
process VARCALL {
    container = 'dancooke/octopus:latest'
    tag "$reference $bamFile"
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/VARCALL"
	debug true
    cpus 1
//    errorStrategy 'ignore'
	
    input:
    path reference
    tuple val(sid), path(bai), path(bamFile)
    path fai
    path bedfile

    output:
    tuple val(sid), path("${sid}.vcf.gz"),      emit:vcf

    script:
    def bed_option = bedfile.getBaseName() == 'dummy' ? "" : "--evaluation-regions ${bedfile}"    // If the base name of bedfile is 'dummy', set bed_option to an empty string

    """    
    octopus \
    -R $reference \
    -I $bamFile \
    --sequence-error-model PCR \
    --forest /opt/octopus/resources/forests/germline.v0.7.4.forest \
    -o ${sid}.vcf.gz \
    --threads ${task.cpus} ${bed_option} \
    -C cancer \
    --very-fast
    """
}