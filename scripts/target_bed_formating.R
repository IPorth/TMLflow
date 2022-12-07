# Script taking input from snakemake workflow and formats bed files for downstream 
# applications Mutect2, CNVkit, Vardict and ONCOCNV

target <- snakemake@input[["bed"]]
#print (target)
Path_ONCOCNV <- snakemake@output[["ONCOCNV"]]
#print(Path_ONCOCNV)
Path_cnvkit <- snakemake@output[["cnvkit"]]
#print(Path_cnvkit)
#Path_vardict <- snakemake@output[["vardict"]]
#print(Path_vardict)


TML_split_bed <- read.table(target, header = FALSE, sep = "\t")

# Target bed for CNVkit and Mutect2

cnvkit_bed <- data.frame(TML_split_bed[,1:3], TML_split_bed[,5])
write.table(cnvkit_bed, Path_cnvkit, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

# Target bed for ONCOCNV

zero <- rep(0, times=length(TML_split_bed[,1]))
Oncocnv_bed <- data.frame(TML_split_bed[,1:4], zero, TML_split_bed[,5])
write.table(Oncocnv_bed, Path_ONCOCNV, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

