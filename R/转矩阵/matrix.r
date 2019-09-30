library(reshape)
t<-read.table("mentel.txt",header=T,sep="\t")
dat2=cast(t,S1~S2,value="NUM")
write.table(dat2,"mentel.xls",sep="\t")
