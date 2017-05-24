#!/usr/bin/env R
args <- commandArgs(trailing=T)
metricsFile <- args[1]
outFile <- args[2]

library(plyr)
#metricsFile <- ("./Total_old_cuff_run/RQCQueryOutput.csv")
##read the information into R
Data <- read.table(metricsFile,header=F,sep=",")
Names <- subset(Data[,1],!duplicated(Data$V1))
Data<- ddply(Data,"V1",numcolwise(sum))
Libraries <- Data[,1]
Data <- Data[,2:7]/3
Data$V1 <- Libraries
Data <- Data[,c(7,1,2,3,4,5,6)]

##set variables for plotting
colors <- c("black","red","blue","orange","green","yellow")
colnames(Data) <- c("name", "Artifacts","JGI.Contam","rRNA","Chloroplast","Mitochondria","E.coli")
#Data <- Data[match(Names,Data[,1]),]
##change the sizes of text and windows for different sizes of input
size <- 1+((96-ifelse(nrow(Data)>96, 96, nrow(Data)))/96)
out.width <- ifelse(nrow(Data)>36, 18, 12)

##write RQC plot to a pdf file
pdf(file=outFile,height=8,width=out.width)
#pdf(file="RQCplot.pdf",height=8,width=out.width)
lapply(2:ncol(Data),function(i){
    if(i==2){
        layout(1:2, heights=c(1, 6))
        par(mar=c(0,0,0,0))
        plot(0, 0, type="n", ann=FALSE, axes=FALSE)
        legend("center", fill=colors, colnames(Data)[-1], cex=1.3, ncol=3)
        par(mar=c(5,5,3,2))
        plot(1:nrow(Data), Data[,i], ylim=c(0,100),
             main="RQC Statistics", xlab="", ylab="",
             xaxt='n', yaxt='n', cex.main=2, cex=size/2,
             pch=19, col=colors[i-1], type='o')
        axis(2, at=seq(0,100, by=10),
             labels=paste(seq(0,100, by=10), "%"), las=2, cex.axis=1.3)
        axis(1, at=seq(1:nrow(Data)),
             labels=Data$name, las=2, cex.axis=1)
        #axis(1,c(48,144,240,336),labels=c("FACS Control","FACS Stressed","C1 Control","C1 Stressed"),cex.axis=1.5)
        #abline(v=c(96.5,192.5,288.5),col="black")
        par(new=TRUE)
        plot(0,0,pch="",ylim=c(0,100), main="", xlab="", ylab="",
             xaxt='n', yaxt='n', abline(h=seq(0,100, by=5), col=rgb(0.5,0.5,1,0.1), lwd=2))
    }
    if(i>2){
        par(new=TRUE)
        plot(1:nrow(Data), Data[,i], ylim=c(0,100),
             main="", xlab="", ylab="",
             xaxt='n', yaxt='n', cex.main=2, cex.lab=1.5, cex=size/2,
             pch=19, col=colors[i-1], type='b')
    }
})   

dev.off()

