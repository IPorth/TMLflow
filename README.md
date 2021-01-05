# TMLflow
Pipeline for TML analysis.
 Uses bam files Ion Torrent S5 sequencer with Tumor Mutation Load Panel.
 Analysis: SNV, CNV, coverage and annotation 

## Step by step
### Install required software
Conda or Miniconda

### Preparing a working directory
Generate a working directory  
`$ mkdir TMLflow`  
`$ cd TMLflow`  
Clone GitHub TMLflow repository into the directory  

### Create and activate the environment
At first time use, the environment hast do be created.  
`$ conda env create --name TMLflow --file environment.yaml`  
Then the environment can be activated by executing:  
`$ conda activate TMLflow`  
For deactiation execute:  
`$ conda deactivate`  

### Required files
- bam files (tumor and normal samples, unmatched)
- target bed file (tab-separated; including chr, start, end, amplicon name and gene name in this order)
- reference genome matching the target bed regions
- adjusted config file 

### Before you start the workflow
To use bwa mem, the reference has to be indexed with the command bwa mem index. An index with samtools faidx is not sufficient.
Execute the following code in the activated TMLflow environment:  
`$ cd refdir`   
`$ bwa mem index ref.fa`  
