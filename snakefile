configfile: "config.yaml"
# Rule 0: includes all files, which should be present at the end of the run.
# Output of the Current last rule is the mapped bam file.
rule all:
    input:
        expand("/mnt/e/TMLflow-merged/mapped/{file}.bam", file=config["files"])

# Rule 1: convertes input bam file to fastq file using picard tools SamToFastq. The step is specific for IONTORRENT data
rule bam_to_fastq:
    input:
        "/mnt/e/TMLflow-merged/data/single/{file}.bam"
    output:
        "/mnt/e/TMLflow-merged/data/fastq/{file}.fastq"
    shell:
        "picard SamToFastq --INPUT {input} --FASTQ {output}"

# Rule 2: mapps fastq generated in rule1 to reference genome hg38
rule bwa_map:
    input:
        "/mnt/e/TMLflow-merged/data/reference/hg38.fa",
        "/mnt/e/TMLflow-merged/data/fastq/{file}.fastq"
    output:
        "/mnt/e/TMLflow-merged/mapped/{file}.bam"
    params:
        rg="@RG\\tID:{file}\\tSM:{file}"
    shell:
        "bwa mem -R '{params.rg}' {input} | samtools view -Sb - > {output}"
