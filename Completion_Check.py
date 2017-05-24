#!/usr/bin/env python

import os
import re
import sys
import smtplib
import pysam

from email.mime.text import MIMEText

Files=os.listdir("./")
Input=sys.argv[1]
InputFile=open(Input,"r")
Last_Header=""
Read_Count=0
Line_Count=0
RNAFilCounter = []
GenomeMapCounter = []
CufflinksCounter = []
TransMapCounter = []
SumCounter = 0
Libraries=[]

for lines in InputFile:
    line=lines.strip("\n")
    line=line.split(" ")
    lib=line[0]
    Libraries.append(lib)
    for File in Files:
        if re.match('^%s$'%lib,File):
            Dir=os.listdir("./%s"%File)
            for subFiles in Dir:
                if re.match("read_stats.txt",subFiles):
                    RS=open("./%s/%s"%(lib,subFiles),"r")
                    for lines in RS:
                        line=lines.strip("\n")
                        line=line.split()
                        Reads=int(line[0])
                        Line_Count += 1
                        if Reads > 0:
                            Read_Count += 1
                    if Read_Count >= 3 and Line_Count >= 3:
                        if lib not in RNAFilCounter:
                            RNAFilCounter.append(lib)
                        else:
                            continue
        if re.match('^%s$'%lib,File):
            Dir=os.listdir("./%s"%File)
            for subFiles in Dir:
                if re.match('%s\.rRNAFiltered.UnpairedRemoved_Read2\.fastq'%lib,subFiles):
                    file_parts=subFiles.split(".")
                    Lib=file_parts[0]
                    Last_Header=file("./%s/%s"%(lib,subFiles),"r").readlines()[-4]
                    Last_Header=Last_Header.strip()
                    Last_Header=Last_Header.split()
                    Last_Header=Last_Header[0]
                    Last_Header=Last_Header[1:]
                if re.match('genome_map',subFiles):
                    Direc=os.listdir("./%s/%s"%(lib,subFiles))
                    for genFiles in Direc:
                        if re.match("accepted_hits.bam",genFiles) or re.match("unmapped.bam",genFiles):
                            mapFile=pysam.AlignmentFile("./%s/%s/%s"%(lib,subFiles,genFiles),"rb")
                            for lines in mapFile:
                                lines=str(lines)
                                line=lines.split("\t")
                                Header=line[0]
                                if Header == str(Last_Header):
                                    if Lib not in GenomeMapCounter:
                                        GenomeMapCounter.append(Lib)
                                    else:
                                        continue
        #if re.match('^%s$'%lib,File):
            
        if re.match('^%s$'%lib,File):
            Dir=os.listdir("./%s"%File)
            for subFiles in Dir:
                if re.match('%s\.rRNAFiltered.UnpairedRemoved_Read2\.fastq'%lib,subFiles):
                    file_parts=subFiles.split(".")
                    Lib=file_parts[0]
                    Last_Header=file("./%s/%s"%(lib,subFiles),"r").readlines()[-4]
                    Last_Header=Last_Header.strip()
                    Last_Header=Last_Header.split()
                    Last_Header=Last_Header[0]
                    Last_Header=Last_Header[1:]
            for subFiles in Dir:
                if re.match('trans_map',subFiles):
                    Direc=os.listdir("./%s/%s"%(lib,subFiles))
                    for genFiles in Direc:
                        if re.match("%s.TransMap.sam"%lib,genFiles):
                            mapFile=open("./%s/%s/%s"%(lib,subFiles,genFiles),"r")
                            for lines in mapFile:
                                lines=lines.strip("\n")
                                line=lines.split("\t")
                                Header=line[0]                      
                                if str(Header) == str(Last_Header):
                                    if Lib not in TransMapCounter:
                                        TransMapCounter.append(Lib)
                                    else:
                                        continue
        if re.match('%s_Sum\.o'%lib,File):
            Sum=open("%s"%File,"r")
            for lines in Sum:
                line=lines.strip("\n")
                if re.match("Done",line):
                    SumCounter += 1

fp=open("./testemail.txt","w+")

if len(Libraries) == len(RNAFilCounter):
   fp.write("All RNAFil Steps Completed"+"\n")
else:
   Missing=list(set(Libraries) - set(RNAFilCounter))
    amount=1
    for items in Missing:
        fp.write("These Libraries did not complete RNAFil:")
        if amount < len(Missing)-1:
            fp.write(items+",")
            amount += 1
        if amount == len(Missing):
            fp.write(items+"\n")
            amount += 1
if len(Libraries) == len(GenomeMapCounter):
   fp.write("All GenomeMap Steps Completed"+"\n")
else:
    Missing=list(set(Libraries) - set(GenomeMapCounter))
    amount=1
    for items in Missing:
        fp.write("These Libraries did not complete GenomeMap:")
        if amount < len(Missing)-1:
            fp.write(items+",")
            amount += 1
        if amount == len(Missing):
            fp.write(items+"\n")
            amount += 1
if len(Libraries) == len(CufflinksCounter):
   fp.write("All Cufflinks Steps Completed"+"\n")
else:
    Missing=list(set(Libraries) - set(CufflinksCounter))
    amount=1
    for items in Missing:
        fp.write("These Libraries did not complete Cufflinks:")
        if amount < len(Missing)-1:
            fp.write(items+",")
            amount += 1
        if amount == len(Missing):
            fp.write(items+"\n")
            amount += 1
if len(Libraries) == len(TransMapCounter):
   fp.write("All TransMap Steps Completed"+"\n")
else:
    Missing=list(set(Libraries) - set(TransMapCounter))
    amount=1
    for items in Missing:
        fp.write("These Libraries did not complete TrnasMap:")
        if amount < len(Missing)-1:
            fp.write(items+",")
            amount += 1
        if amount == len(Missing):
            fp.write(items+"\n")
            amount += 1
if SumCounter == 1:
   fp.write("Sum Step Completed"+"\n")
else:
    print "Sum Step Not Completed"

From=""
To=""

fp.close()
fp=open("./testemail.txt","r")
address=open("./Email.txt","r")
for lines in address:
    line=lines.strip("\n")
    line=line.split("=")
    if re.match("From",lines):
        From=line[1]
    if re.match("To",lines):
        To=line[1]
msg = MIMEText(fp.read())


msg['Subject'] = 'Test for completion check'
msg['From'] = From
msg['To'] = To

s = smtplib.SMTP('localhost')
s.sendmail(From, [To], msg.as_string())
s.quit()




















        #         if re.match('%s\.rRNAFiltered.UnpairedRemoved_Read[0-9]\.fastq'%lib,subFiles):

        # if re.match('%s_GenomeMap\.o'%lib,File):
        #     GenMap=open("%s"%File,"r")
        #     for lines in GenMap:
        #         line=lines.strip("\n")
        #         if re.match("Done",line):
        #             file_parts=File.split("_")
        #             Lib=file_parts[0]
        #             GenomeMapCounter.append(Lib)
        # if re.match('%s_Cuff\.o'%lib,File):
        #     Cuff=open("%s"%File,"r")
        #     for lines in Cuff:
        #         line=lines.strip("\n")
        #         if re.match("Done",line):
        #             file_parts=File.split("_")
        #             Lib=file_parts[0]
        #             CufflinksCounter.append(Lib)
        # if re.match('%s_TransMap\.o'%lib,File):
        #     TransMap=open("%s"%File,"r")
        #     for lines in TransMap:
        #         line=lines.strip("\n")
        #         if re.match("Done",line):
        #             file_parts=File.split("_")
        #             Lib=file_parts[0]
        #             TransMapCounter.append(Lib)
