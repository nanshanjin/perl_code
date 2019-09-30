library(ggplot2)
library(reshape2)
df<-read.table("clean_P1001_CA_DHE06970-18_BHF3JTALXX_L3.QM",sep="\t",head=F)
colnames(df)=c("pos","Q","E")
df$pos<-df$pos+1
middle=(max(df$pos))/2

p <- ggplot(df,aes(x=pos,y=E,ymin=0,ymax=E))+
geom_linerange(colour="#66C2A5",size=0.5)+
xlab("Position along reads") + ylab("% Error rate")+
labs(title="Error rate distribution along reads \n(P1001_CA_DHE06970-18_BHF3JTALXX_L3)")+
geom_vline(xintercept = middle,colour="#377EB8",linetype="dashed")
ggsave(filename="P1001_CA_DHE06970-18_BHF3JTALXX_L3.Error.pdf",plot=p)
ggsave(filename="P1001_CA_DHE06970-18_BHF3JTALXX_L3.Error.png",type="cairo-png",plot=p)
