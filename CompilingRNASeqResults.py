#!/usr/bin/python

import sys

File = sys.argv[1]
Dir=[]
Ref={}
raw={}
cut={}
cutNfil={}

Input = open(File,"r")

for lines in Input:
    line = lines.strip("\n")
    line = line.split(" ")
    Dir.append(line[0])
    Ref[line[-1]]=[]
    Ref[line[-1]].append(line[0])
    
Input.close()

gen={}
gsense={}
ganti={}
gIncSingletons={}
trans={}
tsense={}
tanti={}
tIncSingletons={}


for direc in Dir:
    sunit=direc

    ####### Reading Flagstat Data #############
    GenMapP = open("./%s/genome_map/%s.Primary.GenMap.SORT.flagstat" % (direc,sunit),"r")
    for lines in GenMapP:
        line= lines.strip("\n")
        line= line.split(" ")
        gIncSingletons[direc]=line[0]
        break
    GenMapP.close()
    GenMapPP = open("./%s/genome_map/%s.PP.GenMap.SORT.flagstat" % (direc,sunit),"r")
    for lines in GenMapPP:
        line= lines.strip("\n")
        line= line.split(" ")
        gen[direc]=line[0]
        break
    GenMapPP.close()
    GenMapS = open("./%s/genome_map/%s.PP.GenMap.SORT_sense.flagstat" % (direc,sunit),"r")
    for lines in GenMapS:
        line= lines.strip("\n")
        line= line.split(" ")
        gsense[direc]=line[0]
        break
    GenMapS.close()
    GenMapA = open("./%s/genome_map/%s.PP.GenMap.SORT_antisense.flagstat" % (direc,sunit),"r")
    for lines in GenMapA:
        line= lines.strip("\n")
        line= line.split(" ")
        ganti[direc]=line[0]
        break
    GenMapA.close()
    TransMapP = open("./%s/trans_map/%s.Primary.TransMap.SORT.flagstat" % (direc,sunit),"r")
    for lines in TransMapP:
        line= lines.strip("\n")
        line= line.split(" ")
        tIncSingletons[direc]=line[0]
        break
    TransMapP.close()
    TransMapPP = open("./%s/trans_map/%s.PP.TransMap.SORT.flagstat" % (direc,sunit),"r")
    for lines in TransMapPP:
        line= lines.strip("\n")
        line= line.split(" ")
        trans[direc]=line[0]
        break
    TransMapPP.close()
    TransMapS = open("./%s/trans_map/%s.PP.TransMap.SORT_sense.flagstat" % (direc,sunit),"r")
    for lines in TransMapS:
        line= lines.strip("\n")
        line= line.split(" ")
        tsense[direc]=line[0]
        break
    TransMapS.close()
    TransMapA = open("./%s/trans_map/%s.PP.TransMap.SORT_antisense.flagstat" % (direc,sunit),"r")
    for lines in TransMapA:
        line= lines.strip("\n")
        line= line.split(" ")
        tanti[direc]=line[0]
        break
    TransMapA.close()

    ############# Reading Read Stat Data ###########

    ReadStats = open("./%s/read_stats.txt" % direc,"r")
    Line_Count = 1
    for lines in ReadStats:
        line=lines.strip("\n")
        line=line.split(" ")
        if Line_Count == 1:
            raw[direc]=line[0]
        if Line_Count == 2:
            cut[direc]=line[0]
        if Line_Count == 3:
            cutNfil[direc]=line[0]
        Line_Count+=1
    ReadStats.close()
    
############# Reading Cufflinks Data ##############
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

######## Writing Read Stats ############
RS=open("./Compiled_readstat","w+")
RS.write("Library:,")
i=1
while i <= len(Dir):
    if i < len(Dir):
        RS.write(Dir[i-1]+",")
        i += 1
    elif i == len(Dir):
        RS.write(Dir[i-1])
        i += 1
RS.write("\n")
RS.write("Raw,")
i=1
while i <= len(Dir):
    if i < len(Dir):
        if Dir[i-1] in raw:
            RawCount=float(raw[Dir[i-1]])/4
            RS.write(str(RawCount))
            i += 1
        else: 
            RS.write("")
            i += 1
        RS.write(",")
    elif i == len(Dir):
        if Dir[i-1] in raw:
            RawCount=float(raw[Dir[i-1]])/4
            RS.write(str(RawCount))
            i += 1
        else: 
            RS.write("")
            i += 1
RS.write("\n")
RS.write("AfterTrimming,")
i=1
while i <= len(Dir):
    if i < len(Dir):
        if Dir[i-1] in cut:
            CutCount=float(cut[Dir[i-1]])/4
            RS.write(str(CutCount))
            i += 1
        else:
            RS.write("")
            i += 1
        RS.write(",")
    elif i == len(Dir):
        if Dir[i-1] in cut:
            CutCount=float(cut[Dir[i-1]])/4
            RS.write(str(CutCount))
            i += 1
        else:
            RS.write("")
            i += 1
RS.write("\n")
RS.write("%Trimmed,")
i=1
while i <= len(Dir):
    if i < len(Dir):
        if Dir[i-1] in cut:
            RawCount=float(raw[Dir[i-1]])/4
            CutCount=float(cut[Dir[i-1]])/4
            FilteredPercent=((RawCount-CutCount)/RawCount)*100
            RS.write(str(FilteredPercent))
            i += 1
        else:
            RS.write("")
            i += 1
        RS.write(",")
    elif i == len(Dir):
        if Dir[i-1] in cut:
            RawCount=float(raw[Dir[i-1]])/4
            CutCount=float(cut[Dir[i-1]])/4
            FilteredPercent=((RawCount-CutCount)/RawCount)*100
            RS.write(str(FilteredPercent))
            i += 1
        else:
            RS.write("")
            i += 1
RS.write("\n")
RS.write("After:rRNAFiltering,")
i=1
while i <= len(Dir):
    if i < len(Dir):
        if Dir[i-1] in cutNfil:
            FilCount=float(cutNfil[Dir[i-1]])/4
            RS.write(str(FilCount))
            i += 1
        else:
            RS.write("")
            i += 1
        RS.write(",")
    elif i == len(Dir):
        if Dir[i-1] in cutNfil:
            FilCount=float(cutNfil[Dir[i-1]])/4
            RS.write(str(FilCount))
            i += 1
        else:
            RS.write("")
            i += 1
RS.write("\n")
RS.write("%rRNA,")
i=1
while i <= len(Dir):
    if i < len(Dir):
        if Dir[i-1] in cutNfil:
            RawCount=float(raw[Dir[i-1]])/4
            CutCount=float(cut[Dir[i-1]])/4
            FilCount=float(cutNfil[Dir[i-1]])/4
            FilPercent=((CutCount-FilCount)/RawCount)*100
            RS.write(str(FilPercent))
            i += 1
        else:
            RS.write("")
            i += 1
        RS.write(",")
    elif i == len(Dir):
        if Dir[i-1] in cutNfil:
            RawCount=float(raw[Dir[i-1]])/4
            CutCount=float(cut[Dir[i-1]])/4
            FilCount=float(cutNfil[Dir[i-1]])/4
            FilPercent=((CutCount-FilCount)/RawCount)*100
            RS.write(str(FilPercent))
            i += 1
        else:
            RS.write("")
            i += 1
RS.write("\n")
RS.close()

######## Writing Flagstat #############
Flag=open("./Compiled_flagstat","w+")
Flag.write("LibraryNames:,")
i=1
while i <= len(Dir):
    if i < len(Dir):
        Flag.write(Dir[i-1]+",")
        i += 1
    elif i == len(Dir):
        Flag.write(Dir[i-1])
        i += 1        
Flag.write("\n")
Flag.write("GenMapwithSingletons,")
i=1
while i <= len(Dir):
    if i < len(Dir):
        if Dir[i-1] in gIncSingletons:
            Flag.write(gIncSingletons[Dir[i-1]]+",")
            i += 1
        else:
            Flag.write(",")
            i += 1
    elif i == len(Dir):
        if Dir[i-1] in gIncSingletons:
            Flag.write(gIncSingletons[Dir[i-1]])
            i += 1
        else:
            i += 1        
Flag.write("\n")
Flag.write("%GenMapwithSingletons,")
i=1
while i <= len(Dir):
    if i < len(Dir):
        if Dir[i-1] in gIncSingletons:
            gInc=float(gIncSingletons[Dir[i-1]])
            FilCount=float(cutNfil[Dir[i-1]])/4
            GenMapSing=(gInc/FilCount)*100
            Flag.write(str(GenMapSing)+",")
            i += 1
        else:
            Flag.write(",")
            i += 1
    elif i == len(Dir):
        if Dir[i-1] in gIncSingletons:
            gInc=float(gIncSingletons[Dir[i-1]])
            FilCount=float(cutNfil[Dir[i-1]])/4
            GenMapSing=(gInc/FilCount)*100
            Flag.write(str(GenMapSing))
            i += 1
        else:
            i += 1        
Flag.write("\n")
Flag.write("GenMapPropPaired,")
i=1
while i <= len(Dir):
    if i < len(Dir):
        if Dir[i-1] in gen:
            Flag.write(gen[Dir[i-1]]+",")
            i += 1
        else:
            Flag.write(",")
            i +=1
    elif i == len(Dir):
        if Dir[i-1] in gen:
            Flag.write(gen[Dir[i-1]])
            i +=1
        else:
            i +=1        
Flag.write("\n")
Flag.write("%GenMapPropPaired,")
i=1
while i <= len(Dir):
    if i < len(Dir):
        if Dir[i-1] in gen:
            GenPP=float(gen[Dir[i-1]])
            FilCount=float(cutNfil[Dir[i-1]])/4
            GenMapPP=(GenPP/FilCount)*100
            Flag.write(str(GenMapPP)+",")
            i +=1
        else:
            Flag.write(",")
            i +=1
    elif i == len(Dir):
        if Dir[i-1] in gen:
            GenPP=float(gen[Dir[i-1]])
            FilCount=float(cutNfil[Dir[i-1]])/4
            GenMapPP=(GenPP/FilCount)*100
            Flag.write(str(GenMapPP))
            i +=1
        else:
            i +=1
Flag.write("\n")
Flag.write("GenAntiSensePropPaired,")
i=1
while i <= len(Dir):
    if i < len(Dir):
        if Dir[i-1] in ganti:
            Flag.write(ganti[Dir[i-1]]+",")
            i += 1
        else:
            Flag.write(",")
            i += 1
    elif i == len(Dir):
        if Dir[i-1] in ganti:
            Flag.write(ganti[Dir[i-1]])
            i += 1
        else:
            i += 1        
Flag.write("\n")
Flag.write("%GenAntiSensePropPaired,")
i=1
while i <= len(Dir):
    if i < len(Dir):
        if Dir[i-1] in ganti:
            Ganti=float(ganti[Dir[i-1]])
            Gen=float(gen[Dir[i-1]])
            GantiPer=(Ganti/Gen)*100
            Flag.write(str(GantiPer)+",")
            i += 1
        else:
            Flag.write(",")
            i += 1
    elif i == len(Dir):
        if Dir[i-1] in ganti:
            Ganti=float(ganti[Dir[i-1]])
            Gen=float(gen[Dir[i-1]])
            GantiPer=(Ganti/Gen)*100
            Flag.write(str(GantiPer))
            i += 1
        else:
            i += 1
Flag.write("\n")
Flag.write("TransMapw/tSingletons,")
i=1
while i <= len(Dir):
    if i < len(Dir):
        if Dir[i-1] in tIncSingletons:
            Flag.write(tIncSingletons[Dir[i-1]]+",")
            i +=1
        else:
            Flag.write(",")
            i +=1
    elif i == len(Dir):
        if Dir[i-1] in tIncSingletons:
            Flag.write(tIncSingletons[Dir[i-1]])
            i +=1
        else:
            i +=1        
Flag.write("\n")
Flag.write("%TransMapw/tSingletons,")
i=1
while i <= len(Dir):
    if i < len(Dir):
        if Dir[i-1] in tIncSingletons:
            tInc=float(tIncSingletons[Dir[i-1]])
            FilCount=float(cutNfil[Dir[i-1]])/4
            TransMapPer=(tInc/FilCount)*100
            Flag.write(str(TransMapPer)+",")
            i += 1
        else:
            Flag.write(",")
            i += 1
    elif i == len(Dir):
        if Dir[i-1] in tIncSingletons:
            tInc=float(tIncSingletons[Dir[i-1]])
            FilCount=float(cutNfil[Dir[i-1]])/4
            TransMapPer=(tInc/FilCount)*100
            Flag.write(str(TransMapPer))
            i += 1
        else:
            i += 1        
Flag.write("\n")
Flag.write("TransMapPropPaired,")
i=1
while i <= len(Dir):
    if i < len(Dir):
        if Dir[i-1] in trans:
            Flag.write(trans[Dir[i-1]]+",")
            i +=1
        else:
            Flag.write(",")
            i +=1
    elif i == len(Dir):
        if Dir[i-1] in trans:
            Flag.write(trans[Dir[i-1]])
            i +=1
        else:
            i +=1        
Flag.write("\n")
Flag.write("%TransMapPropPaired,")
i=1
while i <= len(Dir):
    if i < len(Dir):
        if Dir[i-1] in trans:
            TransPP=float(trans[Dir[i-1]])
            FilCount=float(cutNfil[Dir[i-1]])/4
            TransMapPer=(TransPP/FilCount)*100
            Flag.write(str(TransMapPer)+",")
            i += 1
        else:
            Flag.write(",")
            i +=1
    elif i == len(Dir):
        if Dir[i-1] in trans:
            TransPP=float(trans[Dir[i-1]])
            FilCount=float(cutNfil[Dir[i-1]])/4
            TransMapPer=(TransPP/FilCount)*100
            Flag.write(str(TransMapPer))
            i += 1
        else:
            i +=1
Flag.write("\n")
Flag.write("TransSensePropPaired,")
i=1
while i <= len(Dir):
    if i < len(Dir):
        if Dir[i-1] in tsense:
            Flag.write(tsense[Dir[i-1]]+",")
            i += 1
        else:
            Flag.write(",")
            i += 1
    elif i == len(Dir):
        if Dir[i-1] in tsense:
            Flag.write(tsense[Dir[i-1]])
            i += 1
        else:
            i += 1
Flag.write("\n")
Flag.write("%TransSensePropPaired,")
i=1
while i <= len(Dir):
    if i < len(Dir):
        if Dir[i-1] in tsense:
            tSense=float(tsense[Dir[i-1]])
            Trans=float(trans[Dir[i-1]])/4
            tSensePer=(tSense/Trans)*100
            Flag.write(str(tSensePer)+",")
            i += 1
        else:
            Flag.write(",")
            i +=1
    elif i == len(Dir):
        if Dir[i-1] in tsense:
            tSense=float(tsense[Dir[i-1]])
            Trans=float(trans[Dir[i-1]])/4
            tSensePer=(tSense/Trans)*100
            Flag.write(str(tSensePer))
            i += 1
        else:
            i +=1
Flag.write("\n")
Flag.write("TransAntiSensePropPaired,")
i=1
while i <= len(Dir):
    if i < len(Dir):
        if Dir[i-1] in tanti:
            Flag.write(tanti[Dir[i-1]]+",")
            i +=1
        else:
            Flag.write(",")
            i +=1
    elif i == len(Dir):
        if Dir[i-1] in tanti:
            Flag.write(tanti[Dir[i-1]])
            i +=1
        else:
            i +=1
Flag.write("\n")
Flag.write("%TransAntiSensePropPaired,")
i=1
while i <= len(Dir):
    if i < len(Dir):
        if Dir[i-1] in tanti:
            tAnti=float(tanti[Dir[i-1]])
            Trans=float(trans[Dir[i-1]])/4
            tantiPer=(tAnti/Trans)*100
            Flag.write(str(tantiPer)+",")
            i += 1
        else:
            Flag.write(",")
            i +=1
    elif i == len(Dir):
        if Dir[i-1] in tanti:
            tAnti=float(tanti[Dir[i-1]])
            Trans=float(trans[Dir[i-1]])/4
            tantiPer=(tAnti/Trans)*100
            Flag.write(str(tantiPer))
            i += 1
        else:
            i +=1
Flag.write("\n")
Flag.close()
print "Done"
