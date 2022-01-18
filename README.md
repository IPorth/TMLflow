# TMLflow
Pipeline for the analysis of the Oncomine Tumor Mutation Load Assay an Ion AmliSeq Panel.  
It uses bam files from Ion Torrent S5 sequencer with Tumor Mutation Load Panel from tumor ad healthy control FFPE samples.  
The pipeline is designed to run with two replicates, analyses them separately and takes the SNV
found in both replicates with 3% and coverage >100 as output. It can handle more replicates but for the intersection part it will only take replicate 1 and 2. One replicate is not supported yet. 
Analysis: SNV, CNV, coverage, annotation and heterogeneity index.

## Required input files
- bam files (IonTorrent, tumor and normal samples, unmatched, 2 replicates)
- target bed file (tab-separated; including chr, start, end, amplicon name and gene name in this order)
- reference genome matching the target bed regions
- adjusted config file and samples.tsv

## Software and hardware requirements  
Conda or Miniconda
OS: Linux, Windows 10 with WSL2, not tested on Mac
Hardware: Idk, find out!


### Preparations
#### Install required software
Conda or Miniconda  
For information about conda please visit: https://docs.conda.io/en/latest/

#### Preparing a working directory
Generate a working directory  
`mkdir TMLflow`  
`cd TMLflow`  
Clone GitHub TMLflow repository into the directory  

#### Create and activate the environment
At first time, the environment hast to be created.  
`conda env create --name TMLflow --file environment.yaml`  
Then the environment can be activated by executing:  
`conda activate TMLflow`  
For deactiation execute:  
`conda deactivate`  

#### Reference additional files
To use bwa mem, the reference has to be indexed with the command bwa mem index. An index with samtools faidx is not sufficient.
Execute the following code in the activated TMLflow environment:    
`bwa mem index /path/to/ref.fa`  

For Mutect2, a dictionary of your reference has to be generated. This step has to be performed once before you start the workflow for each reference. The reference dictionary has to be in the same directory as the reference itself.   
In the activated environment you can use picard tools to do so:  
`picard CreateSequenceDictionary R=path/to/reference.fasta O=path/to/reference.dict`  

#### Target bed
Bed file requires the columns: chr, start, end, amplicon name, gene name in this particular order.  
Make sure that your target region coordinates are based on the reference use to analyze the data. If not, use liftOver to convert the regions to the right reference (https://genome.ucsc.edu/cgi-bin/hgLiftOver). 

#### Config file and samples.tsv  
Fill in your data information into the samples.tsv. 
Example:  
| bam | sample | condition | rep |
|-----|--------|-----------|-----|
|Sample_A1.bam| Patient_A | Tumor | 1|
|Sample_A2.bam| Patient_A | Tumor | 2|
|Sample_B1.bam| Control_A | Normal | 1|
|Sample_B2.bam| Control_A | Normal | 2|  
  
bam: the exact name of the bam files  
samples: sample name  
condition: either Tumor or Normal, this is case sensitive  
rep: number of replicate  

Insert required location information for the following files and directories into the config file:  
- samples.tsv
- reference, data and output directory
- target bed file
- Annovar directory


#### Download required software and supporting files
##### Annovar
Register on https://annovar.openbioinformatics.org/en/latest/ to download Annovar. Check that the reference genome of supporting databases matches the one used to generate the data.  

##### ONCOCNV
ONCOCNV: search for solution which is not too complicated..


##### Mutect2 SNP database
Download the germline database somatic-hg38_af-only-gnomad.hg38.vcf.gz and its index file from https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38?pli=1.
Store the files in /path/to/TMLflow/support.

##### Final folder structure



The workflow uses different environments. It is recommended to create the environments before starting the workflow.  
`snakemake --use-conda --conda-create-envs-only --cores {NumberOfCores}`  

  
### Output files
