
#Rule1: convertes input bam file to fastq file using picard tools SamToFastq. The step is specific for IONTORRENT data
rule bam_to_fastq:
    input:
        "data/single/{file}.bam"
    output:
        "data/fastq/{file}.fastq"
    shell:
        "picard SamToFastq --INPUT {input} --FASTQ {output}"
