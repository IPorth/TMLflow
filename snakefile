from os.path import join
import pandas as pd

# Reading the config file
configfile: "config.yaml"
# Setting path as global variable
REFDIR=config["Reference"]
DATADIR=config["Data"]
TARGETS=config["targt_bed"]
OUTDIR=config["Output"]
ANNOVAR=config["Annovar"]

#Read sample.tsv
#set_index: indexed the dataframe with the column sample without droping the column out of the dataframe
samples = pd.read_table(config["Sample"], sep="\t", dtype=object).set_index(["sample", "condition","rep"], drop=False)
#Separate conditions from pandas
tumor_only=samples[samples['condition'].str.contains('Tumor')]
tumor_only.columns=['bam', 'sample_tumor', 'condition_tumor',"rep_tumor"]
tumor_only.set_index(["sample_tumor","condition_tumor", "rep_tumor"], drop=False)


control_only=samples[samples['condition'].str.contains('Normal')]
control_only.set_index(["sample","condition", "rep"], drop=False)
# Input functions!
"""
Input function get_files gets all bam file path from pandas and sets the wildcards
"""
def get_files(wildcards):
    return join(DATADIR, samples.loc[(wildcards.sample, wildcards.condition, wildcards.rep), "bam"])

def get_normals(wildcards):
    return control_only.loc(wildcards.sample, wildcards.condition, wildcards.rep)

def get_tumor(wildcards):
    return tumor_only.loc(wildcards.sample, wildcards.condition, wildcards.rep)


# Rule 0: includes all files, which should be present at the end of the run.
# Output of the current last rule of each chapter of the code
rule all:
    input:
        expand(OUTDIR+"/bamstats/{units.sample}_{units.condition}_merge_stats", units=samples.itertuples()),
        expand(OUTDIR+"/annotation/{tumor.sample_tumor}_{tumor.condition_tumor}_1_isec.avinput", tumor=tumor_only.itertuples()),
        expand(OUTDIR+"/CNV/ONCOCNV/{tumor.sample_tumor}_{tumor.condition_tumor}_merge.profile.png", tumor=tumor_only.itertuples(),allow_missing=True)

####################
# Data preparation #
####################

rule bam_to_fastq:
    """convertes input bam file to fastq file using picard tools SamToFastq.
    The step is specific for IONTORRENT data"""
    input:
        get_files
    conda:
       "envs/environment.yaml"
    output:
        OUTDIR+"/fastq/{sample}_{condition}_{rep}.fastq"
    shell:
        "picard SamToFastq --INPUT {input} --FASTQ {output}"

rule bwa_map:
    """mapps fastq to reference genome in REFDIR and saves mapped reads as bam file"""
    input:
        REFDIR,
        OUTDIR+"/fastq/{sample}_{condition}_{rep}.fastq"
    conda:
       "envs/environment.yaml"
    output:
        OUTDIR+"/mapped/{sample}_{condition}_{rep}_realign.bam"
    params:
        rg="@RG\\tID:{sample}_{condition}_{rep}\\tSM:{sample}_{condition}_{rep}"
    shell:
        "bwa mem -R '{params.rg}' {input} | samtools view -Sb - > {output}"

rule samtools_sort:
    """Sorts bam according to chr. Preparation for index and merge"""
    input:
        OUTDIR+"/mapped/{sample}_{condition}_{rep}_realign.bam"
    conda:
       "envs/environment.yaml"
    output:
        OUTDIR+"/sorted/{sample}_{condition}_{rep}.sorted.bam"
    shell:
        "samtools sort -O bam {input} -o {output}"

rule samtools_merge:
    """ Duplicates merged into single file. This step is only required for CNV calling."""
    input:
       #expand(OUTDIR+"/sorted/{units.sample}_{units.condition}_{{rep}}.sorted.bam", units=samples.itertuples(), allow_missing=True)
        OUTDIR+"/sorted/{sample}_{condition}_1.sorted.bam",
        OUTDIR+"/sorted/{sample}_{condition}_2.sorted.bam"
    conda:
       "envs/environment.yaml"
    output:
        OUTDIR+"/merged/{sample}_{condition}_merge.bam"
    shell:
        "samtools merge {output} {input}"

# Rule 6a0: Generation of target bed files
"""Takes input bed file and adjusts it to meet requirements of ONCOCNV and Mutect2 """
rule bed_file_construction:
    input:
        bed=TARGETS
    conda:
       "envs/environment.yaml"
    output:
        cnvkit=OUTDIR+"/target_files/Targets_CNVkit_Mutect.bed",
        ONCOCNV=OUTDIR+"/target_files/Targets_ONCOCNV.bed"
    script:
        "scripts/target_bed_formating.R"

### GATK DepthofCoverage to analyze the files
rule CoverageAnalysis:
    input:
       bam=OUTDIR+"/merged/{sample}_{condition}_merge.bam",
       target=OUTDIR+"/target_files/Targets_CNVkit_Mutect.bed",
       bai=OUTDIR+"/merged/{sample}_{condition}_merge.bam.bai"
    conda:
        "envs/environment.yaml"
    params:
       ref=REFDIR
    output:
       Output1=OUTDIR+"/bamstats/{sample}_{condition}_merge_stats",
       Output2=OUTDIR+"/bamstats/{sample}_{condition}_merge_stats.sample_cumulative_coverage_counts",
       Output3=OUTDIR+"/bamstats/{sample}_{condition}_merge_stats.sample_cumulative_coverage_proportions",
       Output4=OUTDIR+"/bamstats/{sample}_{condition}_merge_stats.sample_interval_statistics",
       Output5=OUTDIR+"/bamstats/{sample}_{condition}_merge_stats.sample_interval_summary",
       Output6=OUTDIR+"/bamstats/{sample}_{condition}_merge_stats.sample_statistics",
       Output7=OUTDIR+"/bamstats/{sample}_{condition}_merge_stats.sample_summary"
    shell:
       "gatk DepthOfCoverage \
       -I {input.bam} \
       -L {input.target} \
       -R {params.ref} \
       -O {output.Output1}"


# Rule 5: Indexing the merged bam files
rule samtools_index_merged:
    input:
        OUTDIR+"/merged/{sample}_{condition}_merge.bam"
    conda:
       "envs/environment.yaml"
    output:
        OUTDIR+"/merged/{sample}_{condition}_merge.bam.bai"
    shell:
        "samtools index {input}"

# Rule 5b: Indexing the single bam files
rule samtools_index_single:
    input:
        OUTDIR+"/sorted/{sample}_{condition}_{rep}.sorted.bam"
    conda:
        "envs/environment.yaml"
    output:
        OUTDIR+"/sorted/{sample}_{condition}_{rep}.sorted.bam.bai"
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
        norm=OUTDIR+"/sorted/{sample}_{condition}_1.sorted.bam",
        bai=OUTDIR+"/sorted/{sample}_{condition}_1.sorted.bam.bai"
    threads: 4
    wildcard_constraints:
        condition= '|'.join([re.escape(x) for x in samples.condition if x == 'Normal']),
    #    rep= '|'.join([re.escape(x) for x in samples.rep if x == 1])
    output:
        OUTDIR+"/normals/{sample}_{condition}_1_mutect2.vcf.gz"
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
        vcf=OUTDIR+"/normals/{sample}_{condition}_1_mutect2.vcf.gz",
        ref=REFDIR
    wildcard_constraints:
        condition= '|'.join([re.escape(x) for x in samples.condition if x == 'Normal']),
    #    rep= '|'.join([re.escape(x) for x in samples.rep if x == 1])
    output:
        OUTDIR+"/normals/{sample}_{condition}_1_mutect2_filtered.vcf.gz"
    shell:
        "gatk FilterMutectCalls \
        -V {input.vcf} \
        -R {input.ref} \
        -O {output}"

# Rule 6b: Generate sample-name-mapped
rule sample_map:
    input:
        sample=expand(OUTDIR+"/normals/{units.sample}_{units.condition}_1_mutect2_filtered.vcf.gz",units=control_only.itertuples(), allow_missing=True )
    wildcard_constraints:
        condition= '|'.join([re.escape(x) for x in samples.condition if x == 'Normal']),
        #rep= '|'.join([re.escape(x) for x in samples.rep if x == 1])
    params:
        name=expand("{units.sample}_{units.condition}_1_mutect2_filtered.vcf.gz", units=control_only.itertuples(), allow_missing=True)
    output:
        OUTDIR+"/normals/sample-name-map.xls"
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
rule mutect2_calling:
    input:
        bam=OUTDIR+"/sorted/{sample}_{condition}_{rep}.sorted.bam",
        bai=OUTDIR+"/sorted/{sample}_{condition}_{rep}.sorted.bam.bai",
        ref=REFDIR,
        target=OUTDIR+"/target_files/Targets_CNVkit_Mutect.bed",
        pon=OUTDIR+"/normals/TML_PoN_single.vcf.gz"
    threads: 4
    wildcard_constraints:
        condition= '|'.join([re.escape(x) for x in samples.condition if x == 'Tumor'])
    conda:
        "envs/environment.yaml"
    output:
        OUTDIR+"/mutect/{sample}_{condition}_{rep}_mutect2.vcf.gz"
    shell:
        """gatk Mutect2 -R {input.ref} -I {input.bam} \
        --intervals {input.target} --native-pair-hmm-threads {threads} \
        --panel-of-normals {input.pon} \
        --max-reads-per-alignment-start 0 -O {output}"""


#Rule 7a: filter Mutect2 calls using gatk FilterMutectCalls
# Annotate variant filters, filter for DP > 250  and AF > 0.1. Remove all mutations with filter status
rule mutect2_filtering:
    input:
        OUTDIR+"/mutect/{sample}_{condition}_{rep}_mutect2.vcf.gz"
    params:
        ref=REFDIR
    wildcard_constraints:
        condition= '|'.join([re.escape(x) for x in samples.condition if x == 'Tumor'])
    conda:
        "envs/filtering.yaml"
    output:
        output1=OUTDIR+"/mutect/filtered/{sample}_{condition}_{rep}_mutect2_filtered.vcf.gz",
        output2=OUTDIR+"/mutect/filtered/{sample}_{condition}_{rep}_mutect2_filtered_selected.vcf.gz",
        output3=OUTDIR+"/mutect/filtered/{sample}_{condition}_{rep}_mutect2_filtered_250_01.vcf",
        output4=OUTDIR+"/mutect/filtered/{sample}_{condition}_{rep}_mutect2_filtered_250_01_cleaned.vcf",
        output5=OUTDIR+"/mutect/filtered/{sample}_{condition}_{rep}_mutect2_filtered_PASS.vcf.gz"
    shell:
        #"scripts/filter_mutect_single.sh {input} {output.output1}  {output.output3}  {output.output5} {params.ref}"
        "scripts/filter_mutect_single.sh {input} {output.output1} {output.output2} {output.output3} {output.output4} {output.output5} {params.ref}"

rule dup_intersection:
    input:
        rep1=OUTDIR+"/mutect/filtered/{sample}_{condition}_1_mutect2_filtered_PASS.vcf.gz",
        rep2=OUTDIR+"/mutect/filtered/{sample}_{condition}_2_mutect2_filtered_PASS.vcf.gz"
    conda:
        "envs/filtering.yaml"
    threads: workflow.cores
    params:
        name="{sample}_{condition}",
        path=OUTDIR+"/intersect"
    output:
        out1=OUTDIR+"/intersect/{sample}_{condition}_1_private.vcf",
        out2=OUTDIR+"/intersect/{sample}_{condition}_2_private.vcf",
        out3=OUTDIR+"/intersect/{sample}_{condition}_1_isec.vcf",
        out4=OUTDIR+"/intersect/{sample}_{condition}_2_isec.vcf"
    shell:
        "scripts/bcf_intersect.sh  {input.rep1} {input.rep2} {params.name} {params.path}"

# Annotation
rule annovar:
    input:
        OUTDIR+"/intersect/{sample}_{condition}_1_isec.vcf"
    params:
        annovar=ANNOVAR,
        file_path=OUTDIR+"/annotation/{sample}_{condition}_1_isec",
        name="{sample}_{condition}_1_isec",
        folder= OUTDIR+"/annotation/"
    output:
        out1=OUTDIR+"/annotation/{sample}_{condition}_1_isec.avinput",
        out2=OUTDIR+"/annotation/{sample}_{condition}_1_isec.hg38_multianno.txt",
        out3=OUTDIR+"/annotation/{sample}_{condition}_1_isec.hg38_multianno.vcf"
    shell:
        "scripts/annovar_2.sh {input} {params.file_path} {params.annovar} {params.name} {params.folder}"


######################
# CNV Calling ONCOCNV#
######################

"""#Prepare variables for ONCOCNV, the part above is identical in the TMLflow snakefile. Combine the workflows when we have decided if single calling is supported
or not. """
ext_T= "_Tumor_merge.bam"
tumor_only_no_index=samples[samples['condition'].str.contains('Tumor')]
Tumor_names=tumor_only_no_index["sample"]
Tumor_names_unique = Tumor_names.unique()
output_tu = ["{}{}".format(i,ext_T) for i in Tumor_names_unique]
Tumor_names_all = ','.join(output_tu)

ext_C= "_Normal_merge.bam"
control_only_no_index=samples[samples['condition'].str.contains('Normal')]
Control_names=control_only_no_index["sample"]
Control_names_unique = Control_names.unique()
output_co = ["{}{}".format(i,ext_C) for i in Control_names_unique]
control_names_all = ','.join(output_co)


rule ONCOCNV_merged:
    input:
        tumor_bam_file=expand(OUTDIR+"/merged/{units.sample_tumor}_{units.condition_tumor}_merge.bam",units=tumor_only.itertuples(), allow_missing=True),
        normal_bam_file=expand(OUTDIR+"/merged/{units.sample}_{units.condition}_merge.bam",units=control_only.itertuples(), allow_missing=True)
    params:
        tumor_bam=Tumor_names_all,
        normal_bam=control_names_all,
        ref_dir=REFDIR,
        out_dir= OUTDIR+"/CNV/ONCOCNV",
        tool_dir= "ONCOCNV-master/ONCOCNV-master/src",
        data_dir= OUTDIR+"/merged",
        target_dir= OUTDIR+"/target_files/Targets_ONCOCNV.bed"
    conda:
        "envs/Oncocnv.yaml"
    output:    
        out1=expand(OUTDIR+"/CNV/ONCOCNV/{units.sample_tumor}_{units.condition_tumor}_merge.profile.png",units=tumor_only.itertuples(), allow_missing=True),
        out2=expand(OUTDIR+"/CNV/ONCOCNV/{units.sample_tumor}_{units.condition_tumor}_merge.profile.txt",units=tumor_only.itertuples(), allow_missing=True),
        out3=expand(OUTDIR+"/CNV/ONCOCNV/{units.sample_tumor}_{units.condition_tumor}_merge.summary.txt",units=tumor_only.itertuples(), allow_missing=True)
    shell:
        "bash ONCOCNV-master/ONCOCNV-master/src/ONCOCNV.sh {params.tumor_bam} {params.normal_bam} {params.ref_dir} {params.out_dir} {params.tool_dir} {params.data_dir} {params.target_dir}"


