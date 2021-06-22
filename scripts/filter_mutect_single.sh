#!/bin/bash
set -euo pipefail

input1=$(realpath $1) # $1 mutect single file 
output1=$(realpath $2) # $3 FilterMutectCalls file 
output2=$(realpath $3) 
output3=$(realpath $4)
output4=$(realpath $5)
output5=$(realpath $6)
ref="$7"

echo "Filter Mutect2 calls"
# Annotation of mutect2 calls with filter reasons
gatk FilterMutectCalls -V ${input1} --max-events-in-region 5 -R ${ref} -O ${output1}
#tabix -f -p vcf ${output1}

#select variants
gatk SelectVariants -R ${ref} -V ${output1} --restrict-alleles-to BIALLELIC -O ${output2}



echo "Applying costum filters AF > 0.03 and DP >100"
# remove calls with less than 3% AF or coverage less than 100
gatk VariantFiltration \
-R ${ref} \
-V ${output2} \
-O ${output3} \
--filter-name "low_coverage" \
--filter-expression "DP < 100" \
--genotype-filter-name "low_AF" \
--genotype-filter-expression "AF < 0.03" 

cat ${output3} | grep -v "low_AF" > ${output4}


echo "Starting vcftools to remove filtered calls"
#vcftools remove-filtered all
vcftools --gzvcf ${output4} --remove-filtered-all --recode --recode-INFO-all --stdout | bgzip -c > ${output5}
tabix -f -p vcf ${output5}

 