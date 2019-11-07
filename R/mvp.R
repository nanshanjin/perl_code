library(MVP)
mydata<-read.table("draw2.txt",sep="\t",header=T)
MVP.Report(pig60K,plot.type="m",LOG10=TRUE,threshold=NULL,col=c("dodgerblue4","deepskyblue"), cex=0.7,chr.den.col=NULL,file="jpg",memo="",dpi=300)
