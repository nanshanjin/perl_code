times<-Sys.time() 

library('getopt')
options(bitmapType='cairo')
spec = matrix(c(
    'infile','f',0,'character',
    'outfile','o',0,'character',
    'lab','l',0,'character'
     ), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
    cat(getopt(spec, usage=TRUE));
    cat("Usage example: \n")
    cat("
    Usage example: 

    Usage:
    --infile sweep out file
    --outfile figure name
    --lab ylab name 
    --help      usage
    \n")
    q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$infile) ) { print_usage(spec) }
if ( is.null(opt$outfile) ) { print_usage(spec) }

mydata<-read.table(opt$infile,header = TRUE,sep="\t")

png(paste(opt$outfile,".png",sep=""),height=900*3,width=1600*3,res=72*3)
ggplot(mydata,aes(x=SNP,y=number,fill=Legend))+geom_bar(stat="identity",width=0.5)
dev.off()
pdf(paste(opt$outfile,".pdf",sep=""),height=900,width=1600)
ggplot(mydata,aes(x=SNP,y=number,fill=Legend))+geom_bar(stat="identity",width=0.5)
dev.off()
escaptime=Sys.time()-times;
print("Done!")
print(escaptime) 
