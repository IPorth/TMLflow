# TMLflow
Pipeline for the analysis of SNV and CNV using data from Oncomine Tumor Mutation Load Assay.  
It uses bam files of tumor and normal FFPE samples from Ion Torrent S5 Prime sequencer.  
The pipeline is designed to run with two replicates, analyses SNVS with Mutect2 tumor-only mode  
separately and takes the SNV found in both replicates.  
Default SNV filtering thresholds are AF>0.1 and DP>250. A panel of normals is utlized in SNV calling.  
CNV calling is performed on merged duplicates with ONCOCNV [https://github.com/BoevaLab/ONCOCNV]


# Flowchart TMLflow
![Flowchart TMLflow](https://github.com/IPorth/TMLflow/blob/main/workflow.png?raw=true)

## Required input files
- bam files (IonTorrent, tumor and normal samples, unmatched, 2 replicates)
- target bed file (tab-separated; including chr, start, end, amplicon name and gene name in this order)
- reference genome matching the target bed regions
- adjusted config file
- samples.tsv

## Software requirements  
Conda or Miniconda installed
For information about conda please visit: https://docs.conda.io/en/latest/  
OS: Ubuntu 18.04.5 , Windows 10 with WSL2 Ubuntu 18.04.6

## Preparations before you start!
### Annovar
Register on https://annovar.openbioinformatics.org/en/latest/ to download Annovar. Check that the reference genome of supporting databases matches the one used to generate the data. Add the annovar directory path to the config file.

### Preparing a working directory
Generate a working directory  
`mkdir TMLflow`  
`cd TMLflow`  
Clone GitHub TMLflow repository into the directory.  

### ONCOCNV
Download ONCOCNV from https://github.com/BoevaLab/ONCOCNV/releases/tag/v6.9.
Place ONCOCNV-master folder into TMLflow directory.  
Open ONCOCNV.sh and exchange the lines below "Set path arguments" until "Run commands" (version 6.9 lines 41-58) with the code in scripts/ONCOCNV_adjust.sh.  
This allows TMLflow to set paths and variables automatically.

### Create and activate the environment
First, the environment has to be created.  
`conda env create --name TMLflow --file environment.yaml`  
Then the environment can be activated by executing:  
`conda activate TMLflow`  
For deactiation execute:  
`conda deactivate`  

### Generate reference additional files
To use bwa mem, the reference has to be indexed with the command bwa mem index. An index with samtools faidx is not sufficient.
Execute the following code in the activated TMLflow environment:    
`bwa mem index /path/to/ref.fa`  

For Mutect2, a dictionary of your reference has to be generated. This step has to be performed once before you start the workflow for each reference. The reference dictionary has to be in the same directory as the reference itself.   
In the activated TMLflow environment you can use picard tools to do so:  
`picard CreateSequenceDictionary R=path/to/ref.fa O=path/to/ref.dict`  

### Target bed
Bed file requires the columns: chr, start, end, amplicon name, gene name in this particular order.  
Make sure that your target region coordinates are based on the reference use to analyze the data. If not, use liftOver to convert the regions to the right reference (https://genome.ucsc.edu/cgi-bin/hgLiftOver).

### Config file and samples.tsv  
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


## Run TMLflow
Make sure the TMLflow environment is activated.   
The workflow uses different environments. It is recommended to create the environments before starting the workflow.  
`snakemake --use-conda --conda-create-envs-only --cores {NumberOfCores}`  
Start the workflow with:  
`snakemake --use-conda --cores {NumberOfCores}`  


## Changing parameters: AF and DP threshold
Open scripts/filter_mutect_single.sh  
DP threshold is set in line 27:  
`--filter-expression "DP < 250`  
Exchange 250 with the desired threshold e.g.  
`--filter-expression "DP < 100`   
Change AF threshold is set in line 29 of the same script in similar way  
Here AF is set to 0.1  
`--genotype-filter-expression "AF < 0.1"`  
Can be changed to e.g. 0.01 or any other threshold:  
`--genotype-filter-expression "AF < 0.01"`  
