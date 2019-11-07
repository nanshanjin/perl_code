
#setwd("/TJPROJ2/HUMAN/Cancer/WES.C101SC16122176.shantoudaxue10lishizhen_cancer.20170105/Mutation/merge_vcf")
dat<-read.table("merged.chr1.snp.xls",head=T,sep="\t",row.names=1)
dat<-data.matrix(dat)
d<-dist(t(dat))
h<-hclust(d,method='ward.D')
png("Samples.cluster.png",width=min(ncol(dat)*50+200,3000),height=600,res=72*2,type='cairo-png')
plot(h,hang=-1,xlab="",ylab="",sub="")
dev.off()
