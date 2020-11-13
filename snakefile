# Reading the config file
configfile: "config.yaml"
# Setting path as global variable
REFDIR=config["Reference"]
DATADIR=config["Data"]


# Rule 0: includes all files, which should be present at the end of the run.
# Output of the Current last rule is the merged bam file.
rule all:
    input:
        expand(DATADIR+"/merged/{sample}_merge.bam.bai", sample=config["samples"]),
        expand(DATADIR+"/normals/{normal}_mutect2.vcf.gz", normal=config["Normals"])


# Rule 1: convertes input bam file to fastq file using picard tools SamToFastq.
# The step is specific for IONTORRENT data
rule bam_to_fastq:
    input:
        DATADIR+"/single/{file}.bam"
    output:
        DATADIR+"/fastq/{file}.fastq"
    shell:
        "picard SamToFastq --INPUT {input} --FASTQ {output}"

# Rule 2: fastq generated in rule 1 mapped to reference genome hg38 in the first step
# Second step saves the mapped reads as bam file and
rule bwa_map:
    input:
        REFDIR,
        DATADIR+"/fastq/{file}.fastq"
    output:
        DATADIR+"/mapped/{file}.bam"
    params:
        rg="@RG\\tID:{file}\\tSM:{file}"
    shell:
        "bwa mem -R '{params.rg}' {input} | samtools view -Sb - > {output}"

# Rule 3:  Sort the bam file from rule 2 according to chr.
# Preparation for index and merge
rule samtools_sort:
    input:
        DATADIR+"/mapped/{file}.bam"
    output:
        DATADIR+"/sorted/{file}.sorted.bam"
    shell:
        "samtools sort -O bam {input} -o {output}"

# Rule 4: All sorted bam files belonging to the same sample are merged into a single file
# This step is not required if there are no replicates
rule samtools_merge:
    input:
        DATADIR+"/sorted/{sample}_1.sorted.bam",
        DATADIR+"/sorted/{sample}_2.sorted.bam"
    output:
        DATADIR+"/merged/{sample}_merge.bam"
    shell:
        "samtools merge {output} {input}"

# Rule 5: Indexing the (merged) bam files
rule samtools_index:
    input:
        DATADIR+"/merged/{sample}_merge.bam"
    output:
        DATADIR+"/merged/{sample}_merge.bam.bai"
    shell:
        "samtools index {input}"

# Rule 6a:
rule mutect2_normal:
    input:
        ref=REFDIR,
        norm=DATADIR+"/merged/{normal}_merge.bam",
        bai=DATADIR+"/merged/{normal}_merge.bam.bai"
    output:
        DATADIR+"/normals/{normal}_mutect2.vcf.gz"
    shell:
        "gatk Mutect2 \
        -R {input.ref} \
        -I {input.norm} \
        -max-mnp-distance 0\
        -O {output}"




# Rule 6: Variant calling with Mutect2 for tumor samples
#rule mutect2_calling:
#    input:
#       bam: DATADIR+"/merged/{tumor}_merge.bam"
#       ref: REFDIR
#       germ:  "af-only-gnomad.vcf.gz"
#       pon:
#
#    output:
#        DATADIR+"/calls/{tumor}_mutect2.vcf"
#    shell:
#          "gatk Mutect2 \
#          -R {input.ref} \
#          -I {input.bam} \
#          -L {}
#          --germline-resource {input.germ} \
#          --panel-of-normals {input.pon}\
#          -O single_sample.vcf.gz"
