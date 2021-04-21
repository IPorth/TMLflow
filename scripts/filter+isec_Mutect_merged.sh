#!/bin/bash
set -euo pipefail

input1=$(realpath $1) # $1 vardict single file 1
output1=$(realpath $2) # $3 FilterMutectCalls file 1
output2=$(realpath $3) # $5 vcftools filtered and zipped file 1
ref=$(realpath $4)


# Filtering mutect
gatk FilterMutectCalls -V ${input1}  -R ${ref} -O ${output1}
tabix -f -p vcf ${output1}

#vcftools remove-filtered all
vcftools --gzvcf ${output1} --remove-filtered-all --recode --recode-INFO-all --stdout | bgzip -c > ${output2}
tabix -f -p vcf ${output2}

