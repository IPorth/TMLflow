configfile: "config.yaml"
# Rule 0: includes all files, which should be present at the end of the run.
# Output of the Current last rule is the merged bam file.
rule all:
    input:
        expand("/mnt/e/TMLflow-merged/data/merged/{sample}_merge.bam", sample=config["samples"])

# Rule 1: convertes input bam file to fastq file using picard tools SamToFastq.
# The step is specific for IONTORRENT data
rule bam_to_fastq:
    input:
        "/mnt/e/TMLflow-merged/data/single/{file}.bam"
    output:
        "/mnt/e/TMLflow-merged/data/fastq/{file}.fastq"
    shell:
        "picard SamToFastq --INPUT {input} --FASTQ {output}"

# Rule 2: fastq generated in rule 1 mapped to reference genome hg38 in the first step
# Second step saves the mapped reads as bam file and
rule bwa_map:
    input:
        "/mnt/e/TMLflow-merged/data/reference/hg38.fa",
        "/mnt/e/TMLflow-merged/data/fastq/{file}.fastq"
    output:
        "/mnt/e/TMLflow-merged/data/mapped/{file}.bam"
    params:
        rg="@RG\\tID:{file}\\tSM:{file}"
    shell:
        "bwa mem -R '{params.rg}' {input} | samtools view -Sb - > {output}"

# Rule 3:  Sort the bam file from rule 2 according to chr.
# Preparation for index and merge
rule samtools_sort:
    input:
        "/mnt/e/TMLflow-merged/data/mapped/{file}.bam"
    output:
        "/mnt/e/TMLflow-merged/data/sorted/{file}.sorted.bam"
    shell:
        "samtools sort -O bam {input} -o {output}"

#Rule 4: All sorted bam files belonging to the same sample are merged into a single file
# This step is not required if there are no replicates
rule samtools_merge:
    input:
        "/mnt/e/TMLflow-merged/data/sorted/{sample}_1.sorted.bam",
        "/mnt/e/TMLflow-merged/data/sorted/{sample}_2.sorted.bam"
    output:
        "/mnt/e/TMLflow-merged/data/merged/{sample}_merge.bam"
    shell:
        "samtools merge {output} {input}"
