#!/usr/bin/python

import json
import sys
import urllib2

Input=sys.argv[1]
Link='https://rqc.jgi-psf.org/api/contam/report/'
libs=[]
InputFile=open(Input,"r")
OutFile=open("./RQCcontam.csv","w+")
for lines in InputFile:
    line=lines.strip("\n")
    line=line.split(" ")
    libs.append(line[0])
InputFile.close()

Lib_list=",".join(libs)
Link=Link+Lib_list

contam_response=urllib2.urlopen(Link)
data=json.load(contam_response)

Libraries=[]
Library=""
artifact=""
contaminants=""
rrna=""
plastid=""
mito=""
ecoli=""

for stat in data['stats']:
    for keys in stat:
        if keys == "library_name":
            library=str(stat[keys])
        if keys == "artifact":
            artifact=str(stat[keys]*100)
        if keys == "contaminants":
            contaminants=str(stat[keys]*100)
        if keys == "rrna":
            rrna=str(stat[keys]*100)
        if keys == "plastid":
            plastid=str(stat[keys]*100)
        if keys == "mitochondrion":
            mito=str(stat[keys]*100)
        if keys == "ecoli_combined":
            ecoli=str(stat[keys]*100)
    OutFile.write(library+","+artifact+","+contaminants+","+rrna+","+plastid+","+mito+","+ecoli+"\n")
        
        

