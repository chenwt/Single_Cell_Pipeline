#!/bin/bash

gff=$1
max_int=$2
ref=$3
name=$4

cd ./$name

echo "Starting: Cufflinks"
/global/homes/k/krgovind/softwares/cufflinks-2.1.1.Linux_x86_64/cufflinks -o cufflink_output -p 16 -I $max_int -G $gff --frag-bias-correct $ref --library-type fr-unstranded ./genome_map/$name.PP.GenMap.SORT.bam
echo "Done"
