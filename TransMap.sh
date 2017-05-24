#!/bin/bash

module load bowtie2
module load samtools
module load picard

name=$1
minins=$2
maxins=$3
t_build=$4

cd ./$name

mkdir -p trans_map
chmod 755 trans_map
echo "Starting: Bowtie-TransMap"
bowtie2 --minins $minins --maxins $maxins  -fr -p 16  -q  -x $t_build -1 ./$name.rRNAFiltered.UnpairedRemoved_Read1.fastq -2 ./$name.rRNAFiltered.UnpairedRemoved_Read2.fastq -S ./trans_map/$name.TransMap.sam
echo "Starting: TransMap Bam Sort"
samtools view -Sb -@ 16 ./trans_map/$name.TransMap.sam > ./trans_map/$name.TransMap.bam
samtools view -Sb -@ 16 ./trans_map/unmapped.sam > ./trans_map/unmapped.bam
samtools sort -@ 16 ./trans_map/$name.TransMap.bam  ./trans_map/$name.TransMap.SORT
samtools index ./trans_map/$name.TransMap.SORT.bam
echo "Starting: TransMap Primary Stats"
samtools view -b -@ 16 -o ./trans_map/$name.Primary.TransMap.SORT.bam -F 256 ./trans_map/$name.TransMap.SORT.bam
samtools flagstat ./trans_map/$name.Primary.TransMap.SORT.bam > ./trans_map/$name.Primary.TransMap.SORT.flagstat
echo "Starting: TransMap Properly Paired Stats"
samtools view -b -@ 16 -o ./trans_map/$name.PP.TransMap.SORT.bam -F 256 -f 2 ./trans_map/$name.TransMap.SORT.bam
echo "Starting: TransMap Strand Splitting"
../split_stranded_rna_seq_reads.py -cf -f ./trans_map/$name.PP.TransMap.SORT.bam
samtools flagstat ./trans_map/$name.PP.TransMap.SORT.bam > ./trans_map/$name.PP.TransMap.SORT.flagstat
samtools flagstat ./trans_map/$name.PP.TransMap.SORT_sense.bam > ./trans_map/$name.PP.TransMap.SORT_sense.flagstat
samtools flagstat ./trans_map/$name.PP.TransMap.SORT_antisense.bam > ./trans_map/$name.PP.TransMap.SORT_antisense.flagstat
samtools index ./trans_map/$name.PP.TransMap.SORT.bam
samtools index ./trans_map/$name.PP.TransMap.SORT_sense.bam
samtools index ./trans_map/$name.PP.TransMap.SORT_antisense.bam
echo "Starting: Picard Insert Size Stats"
picard CollectInsertSizeMetrics H=./trans_map/$name.PP.TransMap.SORT.Hist.pdf I=./trans_map/$name.PP.TransMap.SORT.bam O=./trans_map/$name.PP.TransMap.SORT.Insert.Stats
picard CollectInsertSizeMetrics H=./trans_map/$name.PP.TransMap.SORT_sense.Hist.pdf I=./trans_map/$name.PP.TransMap.SORT_sense.bam O=./trans_map/$name.PP.TransMap.SORT_sense.Insert.Stats
picard CollectInsertSizeMetrics H=./trans_map/$name.PP.TransMap.SORT_antisense.Hist.pdf I=./trans_map/$name.PP.TransMap.SORT_antisense.bam O=./trans_map/$name.PP.TransMap.SORT_antisense.Insert.Stats
echo "Done"
