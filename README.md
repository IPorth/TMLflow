# TMLflow
Pipeline for the analysis of Oncomine Tumor Mutation Load Assay an Ion AmliSeq Panel.  
Uses bam files from Ion Torrent S5 sequencer with Tumor Mutation Load Panel.  
Analysis: SNV, CNV, coverage, annotation and heterogeneity index.

## Required files
- bam files (IonTorrent, tumor and normal samples, unmatched, 2 replicates)
- target bed file (tab-separated; including chr, start, end, amplicon name and gene name in this order)
- reference genome matching the target bed regions
- adjusted config file and metadata.tsv

## Required software  
Conda or Miniconda
OS: Linux, Windows 10 with WSL2, not tested on Mac

## Step by step
### Install required software
Conda or Miniconda

### Preparing a working directory
Generate a working directory  
`mkdir TMLflow`  
`cd TMLflow`  
Clone GitHub TMLflow repository into the directory  

### Create and activate the environment
At first time, the environment hast to be created.  
`conda env create --name TMLflow --file environment.yaml`  
Then the environment can be activated by executing:  
`conda activate TMLflow`  
For deactiation execute:  
`conda deactivate`  

### Before you start the workflow
#### Prepare files
To use bwa mem, the reference has to be indexed with the command bwa mem index. An index with samtools faidx is not sufficient.
Execute the following code in the activated TMLflow environment:  
`cd refdir`   
`bwa mem index ref.fa`  

For Mutect2, a dictionary of your reference has to be generated. This step has to be performed before you start the workflow the first time or when you change your reference.
In the activated environment you can use picard tools to do so:  
`picard CreateSequenceDictionary R=path/reference.fasta O=path/reference.dict`  
  
Fill in your data information into the metadata tsv. Tab separation is required.  
Example:  
| bam | sample | condition | rep |
|-----|--------|-----------|-----|
|Sample_A1.bam| Patient_A | Tumor | 1|
|Sample_A2.bam| Patient_A | Tumor | 2|
|Sample_B1.bam| Control_A | Normal | 1|
|Sample_B2.bam| Control_A | Normal | 2|  
  
bam: the exact name of the input bam files  
samples: sample name  
condition: either Tumor or Normal, the input is case sensitive  
rep: number of replicate  

Insert required information in the config files.  
- Path to metadata tsv
- Path to reference, data and output directory
- Path to target bed file

#### Download required software and 
- Annovar
- ONCOCNV?
- genomAD file for Mutect2 + index file


The workflow uses different environments. It is recommended to create the environments before starting the workflow.  
`snakemake --use-conda --conda-create-envs-only --cores {NumberOfCores}`  

  
  

Somatic Panel of Normals (PoN)  
Please note that once the PoN is created it will not be automaticatlly updated when you add normal samples to the samples.tsv. To include new samples into the PoN, you have to delete the PoN and  the genomicsDB directory and restart the workflow.  
Path to PoN: `output_directory/normals/PoN_single.vcf.gz`  
Path to genomicsDB: `output_directory/pon_db_single`  

