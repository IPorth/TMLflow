#just intersections
#!/bin/bash
set -euo pipefail
input1=$(realpath $1) # path to rep 1 for intersection
input2=$(realpath $2) # path to rep 2 for intersection
name="${3}"
path="${4}" # path where interesected files will be saved

echo "Start intersection of ${name}_1 and ${name}_2! Output path: ${path}"
mkdir -p ${path}
bcftools isec -c some -p ${path} -O v ${input1} ${input2}
cd ${path}
mv 0000.vcf ${name}_1_private.vcf 
mv 0001.vcf ${name}_2_private.vcf
mv 0002.vcf ${name}_1_isec.vcf
mv 0003.vcf ${name}_2_isec.vcf
echo "All done!"