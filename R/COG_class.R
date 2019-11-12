#!/usr/bin/Rscript 
args <- commandArgs(trailingOnly = TRUE)
version <- "v0.1.0, 2017/7/27"

Usage_and_Arguments <- function(author, version){
      cat(paste("
      Author: ", author, "
      Version: ", version, "
      Description: Plot COG class statistics 
      Usage:
           Rscript COG_class.R 
      Example:
           Rscript COG_class.R   all.counts.lwh-4_VS_lwh-1.edgeR_COG_class.xls all.counts.lwh-4_VS_lwh-1.edgeR
"))
}

if ( args[1] == "-h" || args[1] == "-help"){
   Usage_and_Arguments(author, version)
   q()
}

datafile<-args[1]
library(ggplot2)
data<-read.table(datafile,header = T,sep="\t",quote="")
data[,2]<-paste(data[,1],data[,2])
g<-ggplot(data,aes(x=code,y=gene_num,fill=Desciptions))
g<-g+geom_bar(stat="identity")+theme(legend.title=element_blank())+theme(title=element_text(face="bold",size=16,colour="black"),axis.title=element_text(face="bold",size=15,colour="black"),axis.text.x=element_text(face="bold",size=13,colour="black"),axis.text.y=element_text(face="bold",size=15,colour="black"),legend.text=element_text(face="bold",size=13,colour="black"),legend.title=element_text(face="bold",size=13,colour="black"))
g<-g+theme(axis.text.x=element_text(colour="black",size=12,hjust=1))+guides(fill = guide_legend(ncol = 1))
g<-g+labs(x = "", y = "Number of Genes", title = "COG Classification")
ggsave(paste("COG_class.pdf",sep=""),  width=14, height=7)
ggsave(paste("COG_class.png",sep=""),  width=14, height=7)
