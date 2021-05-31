#!/bin/bash
set -euo pipefail

input1=$(realpath $1) # $1 mutect single file 
output1=$(realpath $2) # $3 FilterMutectCalls file 
output2=$(realpath $3) # $5 vcftools filtered and zipped file 
output3=$(realpath $4)
ref=$(realpath $5)

# Annotation of mutect2 calls with filter reasons
gatk FilterMutectCalls -V ${input1}  -R ${ref} -O ${output1}
tabix -f -p vcf ${output1}

# remove calls with less than 3% AF or coverage less than 100
   gatk VariantFiltration \
   -R ${ref} \
   -V ${output1} \
   -O ${output2} \
   --filter-name "Low coverage" \
   --filter-expression "DP < 100" \
   --filter-name "Low AF" \
   --filter-expression "AF < 0.03"

#vcftools remove-filtered all
vcftools --gzvcf ${output2} --remove-filtered-all --recode --recode-INFO-all --stdout | bgzip -c > ${output3}
tabix -f -p vcf ${output3}

