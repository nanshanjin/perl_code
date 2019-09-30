suppressPackageStartupMessages(library("optparse"))
option_list <- list(
        make_option("--infile", action="store", default=NULL, help='The input file.'),
        make_option("--outprefix", action="store", default="Pvalue",help='The outfile prefix.')
	)
#get command line options
opt<-parse_args(OptionParser(usage="%prog [options] file\n", option_list=option_list))
if(is.null(opt$infile))
{ #input file must be given
    cat ("Use  %prog -h for more help info\nThe author: lili@novogene.com")
    quit("no")
}

infile = opt$infile
outfile = paste(opt$outprefix,'.xls',sep="")

data <- read.table(infile,sep= '\t', header=T)
#dat <- as.matrix(data)
##dim 查看维度，dim(dat)表示查看dat维度，100,2表示100行2列
##dim(x) <- c(2,2)表示重新设置的维数2行2列
u <- NULL
for (i in 1:nrow(data)){
    #chr_pos <- data[i,1]
    x <- as.numeric(data[i,5:8])
    dim(x) <- c(2,2)
    a <- fisher.test(x,alternative="two.sided")
    p_orig <- a$p
    or <- as.matrix(a$estimate)[1]
    u   <- rbind(u, data.frame(CHR =data[i,1] ,POS=data[i,2] ,REF=data[i,3] ,ALT=data[i,4] ,Alt_case=data[i,5],Alt_control=data[i,6],Ref_case=data[i,7],Ref_control=data[i,8], Pvalue=p_orig, OR=or))
  }
pvalue <- u[,9]

out <- cbind(u,data.frame(BONF = p.adjust(pvalue,method = "bonferroni" ),HOLM = p.adjust(pvalue,method = "holm" ),FDR_BH =p.adjust(pvalue,method = "BH" ),FDR_BY=p.adjust(pvalue,method = "BY")))

write.table(out,sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE,file=outfile)

