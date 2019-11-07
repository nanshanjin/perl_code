times<-Sys.time()
library('getopt')
library("ComplexHeatmap")
library("circlize")
options(bitmapType='cairo')

spec = matrix(c(
	'm','a',0,'character',
	'o','b',0,'character',
	'help','c', 0, 'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
	
Usage:
	--m different matrix between samples
	--o figure name
	--help		usage
\n")
	q(status=1);
}

if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$m) ) { print_usage(spec) }
if ( is.null(opt$o) ) { print_usage(spec) }



dat_raw <- read.table(opt$m, sep = "\t", header = T,check.names=F)
rownames(dat_raw)<-dat_raw[,1]
dat_raw<-dat_raw[,-1]
dat_raw<-t(dat_raw)

png(paste(opt$o,".png",sep=""),height=900,width=1600)
Heatmap(dat_raw,name = "expression",col = colorRamp2(c(-2, 0, 2),c("green", "white", "red")))
dev.off()
escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
