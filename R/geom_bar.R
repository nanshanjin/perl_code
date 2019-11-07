library("ggplot2")

mydata<-read.delim("A1.txt",sep="\t",header=T)
levels(mydata$sample)
png("A1.png",height=800*3,width=1000*3,res=72*3)#3B-FCF12-CST-3  Hep3B-FCF12-1
ggplot(mydata,aes(Pos,depth,fill=sample))+geom_bar(position="dodge",stat="identity",alpha=0.5)+labs(x = "position", y = "Depth", title = "chr2:136875726-136877225（1.5kb upsteam of CXCR4）",fill="")+theme(title=element_text(face="bold",size=16,colour="black"),axis.title=element_text(face="bold",size=15,colour="black"),axis.text.x=element_text(face="bold",size=13,colour="black"),axis.text.y=element_text(face="bold",size=15,colour="black"),legend.text=element_text(face="bold",size=13,colour="black"),legend.title=element_text(face="bold",size=13,colour="black"))+scale_x_continuous(limits=c(136875726,136877725))+scale_fill_discrete(limits=c("3B-FCF12-CST-3","3B-input5-2-3"))
dev.off()
#scale_fill_manual(values=c("3B-FCF12-CST-3"="Salmon","3B-input5-2-3"="CornflowerBlue"))
pdf("A1.pdf",width=12, height=8)
ggplot(mydata,aes(Pos,depth,fill=sample))+geom_bar(position="dodge",stat="identity",alpha=0.5)+labs(x = "position", y = "Depth", title = "chr2:136875726-136877225（1.5kb upsteam of CXCR4）",fill="")+theme(title=element_text(face="bold",size=16,colour="black"),axis.title=element_text(face="bold",size=15,colour="black"),axis.text.x=element_text(face="bold",size=13,colour="black"),axis.text.y=element_text(face="bold",size=15,colour="black"),legend.text=element_text(face="bold",size=13,colour="black"),legend.title=element_text(face="bold",size=13,colour="black"))+scale_x_continuous(limits=c(136875726,136877725))+scale_fill_discrete(limits=c("3B-FCF12-CST-3","3B-input5-2-3"))
dev.off()
