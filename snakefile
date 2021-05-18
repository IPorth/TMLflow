# Reading the config file
configfile: "config.yaml"
# Setting path as global variable
REFDIR=config["Reference"]
DATADIR=config["Data"]
TARGETS=config["targt_bed"]
OUTDIR=config["Output"]

##JSON controll if all required input files and check path
##Implement conda environments 



# Rule 0: includes all files, which should be present at the end of the run.
# Output of the current last rule of each chapter of the code 
rule all:
    input:
        expand(OUTDIR+"/mutect/MergeEval/{tumor}_2_mutect2_filtered_PASS.vcf.gz", tumor=config["Tumor"]),
        expand(OUTDIR+"/bamstats/{sample}_merge_stats", sample=config["samples"])


####################
# Data preparation #
####################

# Rule 1: convertes input bam file to fastq file using picard tools SamToFastq.
# The step is specific for IONTORRENT data
rule bam_to_fastq:
    input:
        DATADIR+"/{file}.bam"    
    conda:
       "envs/environment.yaml"
    output:
        OUTDIR+"/fastq/{file}.fastq"

    shell:
        "picard SamToFastq --INPUT {input} --FASTQ {output}"

# Rule 2: fastq generated in rule 1 mapped to reference genome hg38 in the first step
# Second step saves the mapped reads as bam file and
rule bwa_map:
    input:
        REFDIR,
        OUTDIR+"/fastq/{file}.fastq"
    conda:
       "envs/environment.yaml"
    output:
        OUTDIR+"/mapped/{file}.bam"
    params:
        rg="@RG\\tID:{file}\\tSM:{file}"
    shell:
        "bwa mem -R '{params.rg}' {input} | samtools view -Sb - > {output}"

# Rule 3:  Sort the bam file from rule 2 according to chr.
# Preparation for index and merge
rule samtools_sort:
    input:
        OUTDIR+"/mapped/{file}.bam"
    conda:
       "envs/environment.yaml"
    output:
        OUTDIR+"/sorted/{file}.sorted.bam"
    shell:
        "samtools sort -O bam {input} -o {output}"

# Rule 4: All sorted bam files belonging to the same sample are merged into a single file
# This step is not required if there are no replicates
rule samtools_merge:
    input:
        OUTDIR+"/sorted/{sample}_1.sorted.bam",
        OUTDIR+"/sorted/{sample}_2.sorted.bam"
    conda:
       "envs/environment.yaml"
    output:
        OUTDIR+"/merged/{sample}_merge.bam"
    shell:
        "samtools merge {output} {input}"

# Rule 6a0: Generation of target bed files  
rule bed_file_construction:
    input:
        bed=TARGETS
    conda:
       "envs/environment.yaml"
    output:
        vardict=OUTDIR+"/target_files/Targets_Vardict.bed",
        cnvkit=OUTDIR+"/target_files/Targets_CNVkit_Mutect.bed",
        ONCOCNV=OUTDIR+"/target_files/Targets_ONCOCNV.bed"
    script:
        "scripts/target_bed_formating.R"

### GATK DepthofCoverage to analyze the files
rule CoverageAnalysis:
    input:
       bam=OUTDIR+"/merged/{sample}_merge.bam",
       target=OUTDIR+"/target_files/Targets_CNVkit_Mutect.bed",
       bai=OUTDIR+"/merged/{sample}_merge.bam.bai"
    conda:
        "envs/environment.yaml"
    params:
       ref=REFDIR
    output:
       Output1=OUTDIR+"/bamstats/{sample}_merge_stats",
       Output2=OUTDIR+"/bamstats/{sample}_merge_stats.sample_cumulative_coverage_counts",
       Output3=OUTDIR+"/bamstats/{sample}_merge_stats.sample_cumulative_coverage_proportions",
       Output4=OUTDIR+"/bamstats/{sample}_merge_stats.sample_interval_statistics",
       Output5=OUTDIR+"/bamstats/{sample}_merge_stats.sample_interval_summary",
       Output6=OUTDIR+"/bamstats/{sample}_merge_stats.sample_statistics",
       Output7=OUTDIR+"/bamstats/{sample}_merge_stats.sample_summary"
    shell:
       "gatk DepthOfCoverage \
       -I {input.bam} \
       -L {input.target} \
       -R {params.ref} \
       -O {output.Output1}"

# Rscript to isolate and formate information

# Rule 5: Indexing the merged bam files
rule samtools_index_merged:
    input:
        OUTDIR+"/merged/{sample}_merge.bam"
    conda:
       "envs/environment.yaml"
    output:
        OUTDIR+"/merged/{sample}_merge.bam.bai"
    shell:
        "samtools index {input}"

# Rule 5b: Indexing the single bam files
rule samtools_index_single:
    input:
        OUTDIR+"/sorted/{file}.sorted.bam"
    conda:
        "envs/environment.yaml"
    output:
        OUTDIR+"/sorted/{file}.sorted.bam.bai"
    shell:
        "samtools index {input}"

###################
#   SNV Calling   #
###################

# Rule 6a1: Generate reference dictionary for use in 6a2
# Picard tools CreateSequenceDictionary --> readme integriert

#Integrate the filtering in the step, if it works do this also for mutect filtering to reduce numer of rules
# Rule 6a2: variant calling for normals, prep for PoN
# No DS
rule mutect2_normal:
    input:
        ref=REFDIR,
        norm=OUTDIR+"/sorted/{normal}_1.sorted.bam",
        bai=OUTDIR+"/sorted/{normal}_1.sorted.bam.bai"
    threads: 4
    output:
        OUTDIR+"/normals/{normal}_1_mutect2.vcf.gz"
    shell:
        "gatk Mutect2 \
        -R {input.ref} \
        -I {input.norm} \
        -max-mnp-distance 0 \
        --max-reads-per-alignment-start 0 \
        --native-pair-hmm-threads {threads} \
        -O {output}"

rule mutect2_normal_filtering:
    input:
        vcf=OUTDIR+"/normals/{normal}_1_mutect2.vcf.gz",
        ref=REFDIR
    output:
        OUTDIR+"/normals/{normal}_1_mutect2_filtered.vcf.gz"
    shell:
        "gatk FilterMutectCalls \
        -V {input.vcf} \
        -R {input.ref} \
        -O {output}"

# Rule 6b: Generate sample-name-mapped
# Wäre cool wenn man diese Regel in Python code umschreiben könnte, damit sie nicht im Flow Diagram auftaucht
rule sample_map:
    input:
        sample=expand(OUTDIR+"/normals/{normal}_1_mutect2_filtered.vcf.gz", normal=config["Normals"])
    output:
        OUTDIR+"/normals/sample-name-map.xls"
    params:
        name=expand("{normal}_1_mutect2_filtered.vcf.gz", normal=config["Normals"]),
    script:
        "scripts/sample-name-map.R"


# Rule 6c: Generation of genomics databas
rule mutect2_GenomicsDB:
    input:
        ref=REFDIR,
        bed=OUTDIR+"/target_files/Targets_CNVkit_Mutect.bed",
        target=OUTDIR+"/target_files/Targets_CNVkit_Mutect.bed",
        normals=OUTDIR+"/normals/sample-name-map.xls"
    threads: workflow.cores
    output:
        directory(OUTDIR+"/pon_db_single")
    shell:
        "gatk GenomicsDBImport -R {input.ref} -L {input.bed} \
        --genomicsdb-workspace-path {output} \
        --validate-sample-name-map TRUE\
        --sample-name-map {input.normals} \
        --intervals {input.target} \
        --merge-input-intervals TRUE "

# Rule 6d: Assemble sommatic panel of normals (PoN)
rule mutect2_PoN_assembyl:
    input:
        ref=REFDIR,
        pon_db=OUTDIR+"/pon_db_single"
    output:
        OUTDIR+"/normals/TML_PoN_single.vcf.gz"
    shell:
        "gatk CreateSomaticPanelOfNormals -R {input.ref} \
        -V gendb://{input.pon_db} -O {output}"

# Muss für Dup1 und Dup2 laufen --> Python script mit allen namen? Wie als wildcard verwenden? Einfach die Regel duplizieren? Meta data sheet?
# Rule 7: SNV calling with Mutect2 for tumor samples
## Dup1
rule mutect2_calling_dup1:
    input:
        bam=OUTDIR+"/sorted/{tumor}_1.sorted.bam",
        ref=REFDIR,
        germ="support/somatic-hg38_af-only-gnomad.hg38.vcf.gz",
        target=OUTDIR+"/target_files/Targets_CNVkit_Mutect.bed",
        bai=OUTDIR+"/sorted/{tumor}_1.sorted.bam.bai",
        pon=OUTDIR+"/normals/TML_PoN_single.vcf.gz"
    threads: 4
    conda:
        "envs/environment.yaml"
    output:
        OUTDIR+"/mutect/{tumor}_1_mutect2.vcf.gz"
    shell:
        """gatk Mutect2 -R {input.ref} -I {input.bam} \
        --intervals {input.target} --native-pair-hmm-threads {threads} \
        --germline-resource {input.germ} \
        --max-reads-per-alignment-start 0 -O {output}"""

# Rule 7a: filter Mutect2 calls using gatk FilterMutectCalls

rule mutect2_filtering_dup1:
    input:
        OUTDIR+"/mutect/{tumor}_1_mutect2.vcf.gz"       
    params:
        ref=REFDIR
    conda:
        "envs/filtering.yaml"
    output:
        output1=OUTDIR+"/mutect/MergeEval/{tumor}_1_mutect2_filtered.vcf.gz",
        output2=OUTDIR+"/mutect/MergeEval/{tumor}_1_mutect2_filtered_PASS.vcf.gz"
    shell:
        "scripts/filter+isec_Mutect_merged.sh {input} {output.output1} {output.output2} {params.ref}"

## Dup2
rule mutect2_calling_dup2:
    input:
        bam=OUTDIR+"/sorted/{tumor}_2.sorted.bam",
        ref=REFDIR,
        germ="support/somatic-hg38_af-only-gnomad.hg38.vcf.gz",
        target=OUTDIR+"/target_files/Targets_CNVkit_Mutect.bed",
        bai=OUTDIR+"/sorted/{tumor}_2.sorted.bam.bai",
        pon=OUTDIR+"/normals/TML_PoN_single.vcf.gz"
    threads: 4
    conda:
        "envs/environment.yaml"
    output:
        OUTDIR+"/mutect/{tumor}_2_mutect2.vcf.gz"
    shell:
        """gatk Mutect2 -R {input.ref} -I {input.bam} \
        --intervals {input.target} --native-pair-hmm-threads {threads} \
        --germline-resource {input.germ} \
        --max-reads-per-alignment-start 0 -O {output}"""

# Rule 7a: filter Mutect2 calls using gatk FilterMutectCalls

rule mutect2_filtering_dup2:
    input:
        OUTDIR+"/mutect/{tumor}_2_mutect2.vcf.gz"       
    params:
        ref=REFDIR
    conda:
        "envs/filtering.yaml"
    output:
        output1=OUTDIR+"/mutect/MergeEval/{tumor}_2_mutect2_filtered.vcf.gz",
        output2=OUTDIR+"/mutect/MergeEval/{tumor}_2_mutect2_filtered_PASS.vcf.gz"
    shell:
        "scripts/filter+isec_Mutect_merged.sh {input} {output.output1} {output.output2} {params.ref}"

# Intersect + analyse VCFS

# Annotation 

# Grafical output SNV Analysis

###############
# CNV Calling #
###############

# DONE separate workflow





