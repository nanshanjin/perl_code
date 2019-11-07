library(ggplot2)
a <- read.table("snp.diff.xls",header=TRUE,sep="\t")
pdf("bar.pdf",width=12, height = 10)
ggplot(a,aes(x=SNP,y=number,fill=Legend))+geom_bar(stat="identity",width=0.5)+theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1))#panel.background = element_blank()È¥µô¿ò
previous_theme <- theme_set(theme_bw)
dev.off()
