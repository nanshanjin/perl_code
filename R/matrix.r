times<-Sys.time()
library("getopt")
library(reshape)
options(bitmapType='cairo')
spec = matrix(c(
    'infile','f',0,'character',
    'outfile','o',0,'character'
    ), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
    cat(getopt(spec, usage=TRUE));
    cat("Usage example: \n")
    cat("
    Usage example: 

    Usage:
    --infile 
    --outfile  out file
    --help      usage
    \n")
    q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$infile) ) { print_usage(spec) }
if ( is.null(opt$outfile) ) { print_usage(spec) }
data<-read.table(opt$infile,header=T,sep="\t")
dat2=cast(data,S1~S2,value="NUM")
write.table(dat2,opt$outfile,sep="\t")
