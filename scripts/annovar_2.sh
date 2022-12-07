#!/bin/bash
#set -euo pipefail

input1=$(realpath $1) # path to rep 1 for intersection
echo $input1
path_out="${2}" # path where interesected files will be saved
echo $path_out
path_annovar="${3}" 
echo $path_annovar
name="${4}" #sample name
echo $name
folder=$(realpath $5) 
echo $folder
#output filenames
#output1=$(realpath $4)
#output2=$(realpath $5)
#output3=$(realpath $6)
#dummy files

linecount=$(grep -v "#" ${input1} | wc -l)
echo $linecount

if [[ $(grep -v -c "#" $input1) -eq 0 ]]
    then
    echo "0 variants in vcf file. Placeholder files generated."
    cd scripts/placeholder
    #kopieren geht nur entweder mit neuem namen im selben directory oder es wird mit gleichem namen in ein neues dir kopiert
    cp placeholder.avinput "$name.avinput" #placeholder1 avinput (empty file)
    mv "$name.avinput" ${folder}
    cp placeholder.txt "$name.hg38_multianno.txt" #placeholder2 txt (column names)
    mv "$name.hg38_multianno.txt" ${folder}
    cp placeholder.vcf "$name.hg38_multianno.vcf" #placeholder3 vcf (annovar vcf header)
    mv "$name.hg38_multianno.vcf" ${folder}
    else
    echo "Number of variants:" $linecount
    cd ${path_annovar}
    perl table_annovar.pl -vcfinput ${input1}  \
    humandb/ -buildver hg38 -out ${path_out} \
    -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation gx,r,f,f,f -nastring . -polish
fi