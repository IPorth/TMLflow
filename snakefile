# Reading the config file
configfile: "config.yaml"
# Setting path as global variable
REFDIR=config["Reference"]
DATADIR=config["Data"]
TARGETS=config["targt_bed"]


# Rule 0: includes all files, which should be present at the end of the run.
# Output of the current last rule of each chapter of the code 
rule all:
    input:
        expand(DATADIR+"/merged/{sample}_merge.bam.bai", sample=config["samples"]),
        expand(DATADIR+"/normals/{normal}_mutect2.vcf.gz", normal=config["Normals"]),
        DATADIR+"/normals/sample-name-map.xls",
        DATADIR+"/pon_db",
        DATADIR+"normals/TML_PoN.vcf.gz",
        expand(DATADIR+"/calls/{tumor}_mutect2.vcf.gz", tumor=config["Tumor"])
        

####################
# Data preparation #
####################


# Rule 1: convertes input bam file to fastq file using picard tools SamToFastq.
# The step is specific for IONTORRENT data
rule bam_to_fastq:
    input:
        DATADIR+"/{file}.bam"
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


###############
# SNV Calling #
###############

# Rule 6a0: Generation of target bed files  
rule bed_file_construction:
    input:
        bed=TARGETS
    output:
        vardict=DATADIR+"/target_files/Targets_Vardict.bed",
        cnvkit=DATADIR+"/target_files/Targets_CNVkit_Mutect.bed",
        ONCOCNV=DATADIR+"/target_files/Targets_ONCOCNV.bed"
    script:
        "scripts/target_bed_formating.R"
# Rule 6a1: Generate reference dictionary for use in 6a2
# Picard tools CreateSequenceDictionary



# Rule 6a2: variant calling for normals, prep for PoN
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
        -max-mnp-distance 0 \
        -O {output}"

# Rule 6b: Generate sample-name-mapped
# Wäre cool wenn man diese Regel in Python code umschreiben könnte, damit sie nicht im Flow Diagram auftaucht
rule sample_map:
    input:
        sample=expand(DATADIR+"/normals/{normal}_mutect2.vcf.gz", normal=config["Normals"])
    output:
        DATADIR+"/normals/sample-name-map.xls"
    params:
        name=expand("{normal}_mutect2.vcf.gz", normal=config["Normals"]),
    script:
        "scripts/sample-name-map.R"


# Rule 6c: Generation of genomics databas
rule mutect2_GenomicsDB :
    input:
        ref=REFDIR,
        bed=DATADIR+"/target_files/Targets_CNVkit_Mutect.bed",
        normals=DATADIR+"/normals/sample-name-map.xls"
    output:
        directory(DATADIR+"/pon_db")
    shell:
        "gatk GenomicsDBImport -R {input.ref} -L {input.bed} \
        --genomicsdb-workspace-path {output} \
        --validate-sample-name-map TRUE\
        --sample-name-map {input.normals} \
        --merge-input-intervals TRUE"



#Rule 6d: Assemble sommatic panel of normals (PoN)
rule mutect2_PoN_assembyl:
    input:
        ref=REFDIR
    params:
        pon_db= "gendb://"+DATADIR+"/pon_db"
    output:
        DATADIR+"/normals/TML_PoN.vcf.gz"
    shell:
        "gatk CreateSomaticPanelOfNormals -R {input.ref} \
        -V {params.pon_db} -O {output}"


#Rule 6: Variant calling with Mutect2 for tumor samples
rule mutect2_calling:
    input:
       bam=DATADIR+"/merged/{tumor}_merge.bam",
       ref=REFDIR,
       germ="support/somatic-hg38_af-only-gnomad.hg38.vcf.gz",
       pon=DATADIR+"/normals/TML_PoN.vcf.gz",
       target=DATADIR+"/target_files/Targets_CNVkit_Mutect.bed"
    output:
        DATADIR+"/calls/{tumor}_mutect2.vcf.gz"
    shell:
          """gatk Mutect2 \
          -R {input.ref} \
          -I {input.bam} \
          -L {input.target} \
          --germline-resource {input.germ} \
          --panel-of-normals {input.pon} \
          --max-reads-per-alignment-start 0 \            
          -O {output} """


#rule vardict:
#    input:
#        ref=REFDIR,
#        target= DATADIR+"/target_files/Targets_Vardict.bed"
#        bam= DATADIR+"/merged/{sample}_merge.bam
#        name="{sample}"
#
#    params:
#        AF_THR= 0.05
#    output:
#
#    shell:
#        "vardict-java -G {input.ref}  -f {params.AF_THR}  -N {input.name}  -b {input.bam} \
#        -c 1 -S 2 -E 3 -g 4 {input.target} | teststrandbias.R | var2vcf_valid.pl -N {input.name} -E -f {params.AF_THR} > {input.name}.vcf"




