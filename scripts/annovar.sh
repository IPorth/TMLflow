#!/bin/bash
set -euo pipefail


input1=$(realpath $1) # path to rep 1 for intersection
path_out="${2}"
path_annovar="${3}" # path where interesected files will be saved

cd ${path_annovar}

perl table_annovar.pl -vcfinput ${input1}  \
humandb/ -buildver hg38 -out ${path_out} \
-remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation gx,r,f,f,f -nastring . -polish