#!/usr/bin/env python

import Bio
from Bio import SeqIO
import sys

fastqFile = sys.argv[1]

fastqBase = fastqFile.replace(".fastq", "")

pair = open(fastqBase+".UnpairedRemoved.fastq","w")
unpair = open(fastqBase+"_unpaired.fq", "w")
handle = open(fastqFile, "rU")
count=0
paired=0
unpaired=0
trimReads={}
for read in SeqIO.parse(handle, "fastq"):
    count+=1
    name = read.description
    baseName = name.split(" ")[0]
    if baseName in trimReads:
        paired +=1
        read1 = trimReads[baseName]
        pair.write(read1.format("fastq"))
        pair.write(read.format("fastq"))
        del trimReads[baseName]
    else:
        trimReads[baseName]=read
    if count % 1000000 ==0:
        print "processed", count/1000000, "M reads"
for name in trimReads:
    unpaired +=1
    read = trimReads[name]
    print >> unpair, read.format("fastq")
print "total paired", paired
print "total unpaired", unpaired
handle.close()
pair.close()
