# octopus
Nextflow-based pipeline for variant calling with octopus


singularity run glebus-sasha-repo-octopus-0.7.4.img \
    -R /home/alexandr/Documents/octopus/results/test/MT.fa \
    -I /home/alexandr/Documents/octopus/results/test/M4-Herk_S41.sorted.bam\
    --config /home/alexandr/Documents/octopus/mitochondria.config \
    --sequence-error-model PCR \
    --forest /opt/octopus/resources/forests/germline.v0.7.4.forest \
    -o 1.vcf.gz \
    --threads 8 