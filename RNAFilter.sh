#!/bin/bash

module load jamo
module load cutadapt
module load duk
module load fastqc
module load bowtie2

################ERCC Spike Ins###########

lib=$1
path=$2
ercc=$3
name=$lib
fastq="$lib.fastq"

cd ./$lib

echo "Starting RNA Filtering:"
jamo link raw_normal library $lib
echo "Starting ZCAT"
zcat *.fastq.gz > $fastq
wc -l $fastq > ./read_stats.txt
echo "Starting Filtering Adapters"
cutadapt -a GATCGGAAGAGCACA -a CTGTCTCTTATA -a AGATCGGAAGAGCGT -O 3 -q 20 -m 32 ./$fastq > ./$name.Cutadapt.fastq
../pe_fastq_splitter.py ./$name.Cutadapt.fastq
wc -l ./$name.Cutadapt.UnpairedRemoved.fastq >> ./read_stats.txt
echo "Starting Filtering rRNA"
duk -k 25 -n ./$name.rRNAFiltered.fastq /global/dna/shared/rqc/ref_databases/qaqc/databases/rRNA.fa ./$name.Cutadapt.UnpairedRemoved.fastq
../pe_fastq_splitter.py ./$name.rRNAFiltered.fastq
wc -l ./$name.rRNAFiltered.UnpairedRemoved.fastq  >> ./read_stats.txt
echo "Starting Separating Reads"
../read_seperator.py ./$name.rRNAFiltered.UnpairedRemoved.fastq
fastqc --noextract -t 3 ./$fastq ./$name.Cutadapt.UnpairedRemoved.fastq ./$name.rRNAFiltered.UnpairedRemoved.fastq
rm -rf ./$name.Cutadapt.fastq ./$name.Cutadapt.UnpairedRemoved.fastq ./$name.rRNAFiltered.fastq
echo "Starting Searching for ERCC Spike-Ins"
bowtie2 --maxins 2000 -fr -p 16 -q -x $ercc -1 "./$lib.rRNAFiltered.UnpairedRemoved_Read1.fastq" -2 "./$lib.rRNAFiltered.UnpairedRemoved_Read2.fastq" -S "./$lib.Spikein.sam"
cat ./$lib.Spikein.sam | grep -v "^@" | cut -f 3 | grep -v "*" | sort | uniq -c > ./$lib.spikeins.bow.txt
cat ./$lib.spikeins.bow.txt | wc -l >> SpikeinNumbersCounts.bow.txt
echo "Done"
