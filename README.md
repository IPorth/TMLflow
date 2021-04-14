# TMLflow
Pipeline for the analysis of Oncomine Tumor Mutation Load Assay an Ion AmliSeq Panel.  
Uses bam files from Ion Torrent S5 sequencer with Tumor Mutation Load Panel.  
Analysis: SNV, CNV, coverage, annotation and heterogeneity index.   

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

The workflow uses different environments. It is recommended to create the environments before starting the workflow.  
`snakemake --use-conda --conda-create-envs-only --cores {NumberOfCores}`

### Required files
- bam files (tumor and normal samples, unmatched)
- target bed file (tab-separated; including chr, start, end, amplicon name and gene name in this order)
- reference genome matching the target bed regions
- adjusted config file 

### Before you start the workflow
#### Prepare files
To use bwa mem, the reference has to be indexed with the command bwa mem index. An index with samtools faidx is not sufficient.
Execute the following code in the activated TMLflow environment:  
`cd refdir`   
`bwa mem index ref.fa`  

For Mutect2, a dictionary of your reference has to be generated. This step has to be performed before you start the workflow the first time or when you change your reference.
In the activated environment you can use picard tools to do so:  
`picard CreateSequenceDictionary R=path/reference.fasta O=path/reference.dict`  

Insert required information in the config files.  
- Normal and tumor sample name 
- Path to reference, data and output directory
- Path to target bed file

#### Download required software and 
- Annovar
- ONCOCNV?
- genomAD file for Mutect2 + index file

