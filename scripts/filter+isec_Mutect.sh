#!/bin/bash
set -euo pipefail

input1=$(realpath $1) # $1 vardict single file 1
input2=$(realpath $2) # $2 vardict single file 2
output1=$(realpath $3) # $3 FilterMutectCalls file 1
output2=$(realpath $4) # $4 FilterMutectCalls file 2
output3=$(realpath $5) # $5 vcftools filtered and zipped file 1
output4=$(realpath $6) # $6 vcftools filtered and zipped file 2
name="$7" 
path="${8}"
ref="$9"

echo "Starting FilterMutectcalls"
# Filtering mutect
gatk FilterMutectCalls -V ${input1}  -R ${ref} -O ${output1}
tabix -f -p vcf ${output1}
gatk FilterMutectCalls -V ${input2}  -R ${ref} -O ${output2}
tabix -f -p vcf ${output2}

echo "Remove all mutations which do not pass filters"
#vcftools remove-filtered all
vcftools --gzvcf ${output1} --remove-filtered-all --recode --recode-INFO-all --stdout | bgzip -c > ${output3}
tabix -f -p vcf ${output3}
vcftools --gzvcf ${output2} --remove-filtered-all --recode --recode-INFO-all --stdout | bgzip -c > ${output4}
tabix -f -p vcf ${output4}

echo "Intersection of ${output3} and ${output4}"
#intersection
echo ${path}
mkdir -p ${path}
bcftools isec -p ${path} -O v ${output3} ${output4}
cd ${path}
mv 0000.vcf ${name}_1_private.vcf 
mv 0001.vcf ${name}_2_private.vcf
mv 0002.vcf ${name}_1_isec.vcf
mv 0003.vcf ${name}_2_isec.vcf
echo "All done!"