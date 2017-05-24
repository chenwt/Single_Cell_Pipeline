#!/usr/bin/python

import sys

File = sys.argv[1]
Dir=[]
Ref={}

Input = open(File,"r")

for lines in Input:
    line = lines.strip("\n")
    line = line.split(" ")
    Dir.append(line[0])
    Ref[line[-1]]=[]
    Ref[line[-1]].append(line[0])
    
Input.close()

gene={}
locus={}
length={}
cov={}
fpkm={}


for key in Ref:
    cuffdirs=Ref[key]
    for direc in cuffdirs:
        Cuff=open("./%s/cufflink_output/isoforms.fpkm_tracking" % direc,"r")
        gene[direc]={}
        locus[direc]={}
        length[direc]={}
        cov[direc]={}
        fpkm[direc]={}
        for lines in Cuff:
            line = lines.strip("\n")
            line = line.split("\t")
            gene[direc][line[0]]=line[3]
            locus[direc][line[0]]=line[6]
            length[direc][line[0]]=line[7]
            cov[direc][line[0]]=line[8]
            fpkm[direc][line[0]]=line[9]
        Cuff.close()
    Out=open("./Compiled_cufflinks.%s" % key,"w+")
    Out.write("Tracking_Id" + "\t" + "Gene_Id" + "\t" + "Locus" + "\t" + "Length" + "\t")
    i=1
    while i <= len(cuffdirs):
        if i < len(cuffdirs):
            Out.write(cuffdirs[i-1] + ".cov" + "\t" + cuffdirs[i-1] + ".fpkm" + "\t")
            i += 1
        elif i == len(cuffdirs):
            Out.write(cuffdirs[i-1] + ".cov" + "\t" + cuffdirs[i-1] + ".fpkm" + "\n")
            i += 1
    for keys in gene[cuffdirs[0]]:
        Out.write(keys + "\t" + gene[cuffdirs[0]][keys] + "\t" + locus[cuffdirs[0]][keys] + "\t" + length[cuffdirs[0]][keys] + "\t")
        i=1
        while i <= len(cuffdirs):
            if i < len(cuffdirs):
                if keys in cov[cuffdirs[i-1]]:
                    Out.write(cov[cuffdirs[i-1]][keys] + "\t")
                else:
                    Out.write("NA" + "\t")
                if keys in fpkm[cuffdirs[i-1]]:
                    Out.write(fpkm[cuffdirs[i-1]][keys] + "\t")
                else:
                    Out.write("NA" + "\t")
                i += 1
            elif i == len(cuffdirs):
                if keys in cov[cuffdirs[i-1]]:
                    Out.write(cov[cuffdirs[i-1]][keys] + "\t")
                else:
                    Out.write("NA" + "\t")
                if keys in fpkm[cuffdirs[i-1]]:
                    Out.write(fpkm[cuffdirs[i-1]][keys] + "\n")
                else:
                    Out.write("NA" + "\n")
                i += 1

    Out.close()
