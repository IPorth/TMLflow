#just intersections
#!/bin/bash
set -euo pipefail
input1=$(realpath $1) # path to file 1 for intersection
input2=$(realpath $2) # path to file 2 for intersection
file1="${3}" # name file 1
file2="${4}" # name file 2
path="${5}" # path where interesected files will be saved

echo ${path}
mkdir -p ${path}
bcftools isec -p ${path} -O v ${input1} ${input2}
cd ${path}
mv 0000.vcf ${file1}_private.vcf 
mv 0001.vcf ${file2}_private.vcf
mv 0002.vcf ${file1}_isec.vcf
mv 0003.vcf ${file2}_isec.vcf
echo "All done!"