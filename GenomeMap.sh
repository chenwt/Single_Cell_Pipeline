#!/bin/bash

module load tophat
module load samtools
module load picard
module load bowtie2

lib=$1
stddev=$2
min_int=$3
max_int=$4
gff=$5
genin=$6
name=$7

cd ./$name

mkdir -p genome_map
#chmod -R 755 genome_map
chmod 755 genome_map
echo "Starting: TophatGenMap"
tophat -o genome_map -r 150 --mate-std-dev $stddev -i $min_int -I $max_int -p 16 -G $gff --library-type fr-unstranded $genin ./$name.rRNAFiltered.UnpairedRemoved_Read1.fastq ./$name.rRNAFiltered.UnpairedRemoved_Read2.fastq
echo "Starting: Bam Filtering and Flagstats"
ln -s accepted_hits.bam ./genome_map/$name.GenMap.bam
samtools sort -@ 16 ./genome_map/$name.GenMap.bam ./genome_map/$name.GenMap.SORT
samtools index ./genome_map/$name.GenMap.SORT.bam
echo "Starting: Primary Alignment Stats"
samtools view -b -@ 16 -o ./genome_map/$name.Primary.GenMap.SORT.bam -F 256 ./genome_map/$name.GenMap.SORT.bam
samtools flagstat ./genome_map/$name.Primary.GenMap.SORT.bam > ./genome_map/$name.Primary.GenMap.SORT.flagstat
echo "Starting: Properly Paired Stats"
samtools view -b -@ 16 -o ./genome_map/$name.PP.GenMap.SORT.bam -F 256 -f 2 ./genome_map/$name.GenMap.SORT.bam
samtools flagstat ./genome_map/$name.PP.GenMap.SORT.bam > ./genome_map/$name.PP.GenMap.SORT.flagstat
echo "Starting: Splitting Strands"
../split_stranded_rna_seq_reads.py -cf -f ./genome_map/$name.PP.GenMap.SORT.bam
samtools flagstat ./genome_map/$name.PP.GenMap.SORT_sense.bam > ./genome_map/$name.PP.GenMap.SORT_sense.flagstat
samtools flagstat ./genome_map/$name.PP.GenMap.SORT_antisense.bam > ./genome_map/$name.PP.GenMap.SORT_antisense.flagstat
samtools index ./genome_map/$name.PP.GenMap.SORT.bam
samtools index ./genome_map/$name.PP.GenMap.SORT_sense.bam
samtools index ./genome_map/$name.PP.GenMap.SORT_antisense.bam
echo "Starting: Picard Insert Size"
picard CollectInsertSizeMetrics H=./genome_map/$name.PP.GenMap.SORT.Hist.pdf I=./genome_map/$name.PP.GenMap.SORT.bam O=./genome_map/$name.PP.GenMap.SORT.Insert.Stats

picard CollectInsertSizeMetrics H=./genome_map/$name.PP.GenMap.SORT_sense.Hist.pdf I=./genome_map/$name.PP.GenMap.SORT_sense.bam O=./genome_map/$name.PP.GenMap.SORT_sense.Insert.Stats

picard CollectInsertSizeMetrics H=./genome_map/$name.PP.GenMap.SORT_antisense.Hist.pdf I=./genome_map/$name.PP.GenMap.SORT_antisense.bam O=./genome_map/$name.PP.GenMap.SORT_antisense.Insert.Stats

rm -rf ./genome_map/$name.GenMap.bam ./genome_map/$name.GenMap.bam.bai ./genome_map/$name.Primary.GenMap.SORT.bam ./genome_map/$name.PP.GenMap.SORT_sense.bam ./genome_map/$name.PP.GenMap.SORT_sense.bam.bai ./genome_map/$name.PP.GenMap.SORT_antisense.bam ./genome_map/$name.PP.GenMap.SORT_antisense.bam.bai
echo "Done"
