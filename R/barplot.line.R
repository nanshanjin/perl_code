times<-Sys.time()
library("getopt")
library("ggplot2")
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
data <- read.table(opt$infile,header=T,sep="\t")
png(paste(opt$outfile,".png",sep=""),height=800*3,width=800*3,res=72*3)
ggplot(data,aes(x=distance,y=..density..))+geom_histogram(fill = "lightblue", colour = "black")+stat_density(geom='line',position='identity',size=1.2)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
pdf(paste(opt$outfile,".pdf",sep=""),height=800,width=800)
ggplot(data,aes(x=distance,y=..density..))+geom_histogram(fill = "lightblue", colour = "black")+stat_density(geom='line',position='identity',size=1.2)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
