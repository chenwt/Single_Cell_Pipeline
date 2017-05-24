#!/usr/bin/python

import sys

File = sys.argv[1]

forward=File
forward=forward.replace(".fastq","_Read1.fastq")

reverse=File
reverse=reverse.replace(".fastq","_Read2.fastq")

Input = open(File,"r")
FOutput = open(forward,"w+")
ROutput = open(reverse,"w+")

Line_Count = 1

for lines in Input:
    line = lines.strip("\n")
    if Line_Count <= 4:
        FOutput.write(line + "\n")
        Line_Count += 1
    elif Line_Count < 8:
        ROutput.write(line + "\n")
        Line_Count += 1
    elif Line_Count :
        ROutput.write(line + "\n")
        Line_Count = 1

