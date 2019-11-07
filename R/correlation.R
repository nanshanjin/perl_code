times<-Sys.time()
library("RColorBrewer",lib.loc="/home/Group/rna/rna-chip/soft/R/")
library("pheatmap",lib.loc="/home/Group/rna/rna-chip/soft/R/")
library("Heatplus",lib.loc="/home/Group/rna/rna-chip/soft/R/")
library('getopt')
options(bitmapType='cairo')
spec = matrix(c(
    'infile','f',0,'character',
    'outfile','o',0,'character'
     ), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
    cat(getopt(spec, usage=TRUE));
    cat("
    Usage:
    --infile :input nordata file
    --outfile :out file prefix
    --help      usage
    \n")
    q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$infile) ) { print_usage(spec) }
if ( is.null(opt$outfile) ) { print_usage(spec) }
mydata<-read.delim(opt$infile,header = TRUE,sep="\t",row.names=1)
####方法一
result <- NULL
for (i in 1:nrow(mydata)){
    data<-as.numeric(mydata[i,])
    a<-(data-mean(data))/sd(data)
    result<-rbind(result,a)
    }
colnames(result)<-colnames(mydata)
rownames(result)<-rownames(mydata)
#z-score 方法二
#mydata.zscore <- apply(mydata, 1, scale)
#result<-t(mydata.zscore)
write.table(result,paste(opt$outfile,"zscore.txt",sep="_"),quote = F,sep = "\t",col.names=NA)
outpdf<-paste(opt$outfile,".pdf",sep="")
myplot <- pheatmap(result,cluster_rows=T,cluster_cols=T, 
    fontsize_number = 8,
    border_color=NA,
        filename=outpdf,
        )
outpng<-paste(opt$outfile,".png",sep="")
myplot <- pheatmap(result,cluster_rows=T,cluster_cols=T, 
fontsize_number = 8,
border_color=NA,
filename=outpng,
)
escaptime=Sys.time()-times;
print("Done!")
print(escaptime) 
