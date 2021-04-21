#!/bin/bash
set -euo pipefail

input1=$(realpath $1) # $1 vardict single file merged
output1=$(realpath $2) # $2 vcftools filtered and zipped file merged
output2=$(realpath $3) # $3 vcffilter filtered file 2 0.05
output3=$(realpath $4) # $4 vcffilter filtered file 2 0.1
name="$5" 




# Filtering vardict 
#file 1
vcftools --vcf ${input1} --remove-filtered-all --recode --recode-INFO-all --stdout | bgzip -c > ${output1}
tabix -f -p vcf ${output1}
vcffilter -f "AF >0.05" ${output1}| bgzip -c > ${output2}
tabix -f -p vcf ${output2}
vcffilter -f "AF >0.1" ${output1} | bgzip -c > ${output3}
tabix -f -p vcf ${output3}
vcffilter -f "AF >0.05" -f "AF <0.1" ${output1} > ${name}_vardict_5-10.vcf

