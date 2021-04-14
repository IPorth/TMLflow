#!/usr/bin/env sh
input1=$(realpath $1) # $1 vardict single file 1
input2=$(realpath $2) # $2 vardict single file 2
output1=$(realpath $3) # $3 vcftools filtered and zipped file 1
output2=$(realpath $4) # $4 vcftools filtered and zipped file 2
output3=$(realpath $5) # $5 vcffilter filtered file 1 0.05
output4=$(realpath $6) # $6 vcffilter filtered file 2 0.05
output5=$(realpath $7) # $7 vcffilter filtered file 1 0.1
output6=$(realpath $8) # $8 vcffilter filtered file 2 0.1
name="$9" 
path=$(realpath ${10})



# Filtering vardict 
#file 1
vcftools --gzvcf ${input1} --remove-filtered-all --recode --recode-INFO-all --stdout | bgzip -c > ${output1}
tabix -p vcf ${output1}
vcffilter -f "AF >0.05" ${output1}| bgzip -c > ${output3}
tabix -p vcf ${output3}
vcffilter -f "AF >0.1" ${output1}| bgzip -c > ${output5}
tabix -p vcf ${output5}
#file 2
vcftools --gzvcf ${input2} --remove-filtered-all --recode --recode-INFO-all --stdout | bgzip -c > ${output2}
tabix -p vcf ${output2}
vcffilter -f "AF >0.05" ${output2}| bgzip -c  > ${output4}
tabix -p vcf ${output4}
vcffilter -f "AF >0.1" ${output2}| bgzip -c > ${output6}
tabix -p vcf ${output6}
#intersection

echo ${path}
bcftools isec -p ${path} -O v ${output1} ${output2}
cd ${path}
mv 0000.vcf ${name}_1_private.vcf 
mv 0001.vcf ${name}_2_private.vcf
mv 0002.vcf ${name}_1_isec.vcf
mv 0003.vcf ${name}_2_isec.vcf

bcftools isec -p ${path} -O v ${output3} ${output4}
cd ${path}
mv 0000.vcf ${name}_1_5_private.vcf 
mv 0001.vcf ${name}_2_5_private.vcf
mv 0002.vcf ${name}_1_5_isec.vcf
mv 0003.vcf ${name}_2_5_isec.vcf

bcftools isec -p ${path} -O v ${output6} ${output5}
cd ${path}
mv 0000.vcf ${name}_1_10_private.vcf 
mv 0001.vcf ${name}_2_10_private.vcf
mv 0002.vcf ${name}_1_10_isec.vcf
mv 0003.vcf ${name}_2_10_isec.vcf