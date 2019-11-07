times<-Sys.time() 
library(ggplot2)
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
    --infile   file
    --outfile figure name
    --lab ylab name 
    --help      usage
    \n")
    q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$infile) ) { print_usage(spec) }
if ( is.null(opt$outfile) ) { print_usage(spec) }

a<-read.table(opt$infile,header = F,sep="\t",stringsAsFactors=F)
g<-sprintf("%s%s","Chr",unique(a[,1]))
a[,7]<-a[,3]-a[,2]
for(i in 1:length(a[,1])){
  if(a[i,4] == "478"){
    a[i,8]=a[i,1]-0.2
  }else{
    a[i,8]=a[i,1]+0.2
  }
} 
colors<-c("red","blue","green")
p<-ggplot(a,aes(x=a[,8],y=a[,7],fill=a[,5],group=a[,4]))+geom_bar(stat = "identity",width = 0.5)+    
  geom_text(aes(label = a[,6]), size = 3, hjust = 0, vjust = 0,position = "stack",angle=45)+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     legend.position = c(0.9,0.9),legend.title=element_blank())+
  scale_fill_manual(values = colors)+
  scale_x_continuous(breaks = c(1:10),labels = g)+labs(x="",y="")+
  scale_y_continuous(labels = c(0,"100MB","200MB","300MB"))
 
pdf(paste(opt$outfile,".pdf",sep=""),height=800,width=800)
print(p)
dev.off()
png(paste(opt$outfile,".png",sep=""),height=800*3,width=800*3,res=72*3)
print(p)
dev.off()

escaptime=Sys.time()-times;
print("Done!")
print(escaptime) 
