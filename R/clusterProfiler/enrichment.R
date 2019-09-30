times<-Sys.time() 

library('getopt')
options(bitmapType='cairo')
spec = matrix(c(
    'input','f',0,'character',
    'outdir','o',0,'character'
     ), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
    cat(getopt(spec, usage=TRUE));
    cat("Usage example: \n")
    cat("
    Usage example: 

    Usage:
    --input input file
    --outdir out dir
    --help      usage
    \n")
    q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$input) ) { print_usage(spec) }
if ( is.null(opt$outdir) ) { print_usage(spec) }

library(clusterProfiler)
library(org.Hs.eg.db)
#bitr fromType:ACCNUM, ALIAS, ENSEMBL, ENSEMBLPROT, ENSEMBLTRANS, ENTREZID, ENZYME, EVIDENCE, EVIDENCEALL, GENENAME, GO, GOALL, IPI, MAP, OMIM, ONTOLOGY, ONTOLOGYALL, PATH, PFAM, PMID, PROSITE, REFSEQ, SYMBOL, UCSCKG, UNIGENE, UNIPROT

#eg <- bitr("ENSG00000282591", fromType="GENENAME", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
mydata <- read.delim(opt$input,sep="\t",header=T)
eg <- bitr(mydata$SYMBOL,fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")

genelist <- eg$ENTREZID
go <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod ='bonferroni',pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.2,keyType = 'ENTREZID')

kegg <- enrichKEGG(genelist, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod ='bonferroni',
                 minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)
setwd(opt$outdir)
write.table(go,file="go_enrichment.txt",sep="\t",quote=F,row.names=F)
write.table(kegg,file="kegg_enrichment.txt",sep="\t",quote=F,row.names=F)
png("go_bar.png",height=800*3,width=800*3,res=72*3)
barplot(go,showCategory=20,drop=T)
dev.off()
pdf("go_bar.pdf",height=800,width=800)
barplot(go,showCategory=20,drop=T)
dev.off()

png("go_dot.png",height=800*3,width=800*3,res=72*3)
dotplot(go,showCategory=20)
dev.off()
pdf("go_dot.pdf",height=800,width=800)
dotplot(go,showCategory=20)
dev.off()

png("kegg_bar.png",height=800*3,width=800*3,res=72*3)
barplot(go,showCategory=20,drop=T)
dev.off()
pdf("kegg_bar.pdf",height=800,width=800)
barplot(go,showCategory=20,drop=T)
dev.off()

png("kegg_dot.png",height=800*3,width=800*3,res=72*3)
dotplot(kegg, showCategory=30)
dev.off()
pdf("kegg_dot.pdf",height=800,width=800)
dotplot(kegg, showCategory=30)
dev.off()

escaptime=Sys.time()-times;
print("Done!")
print(escaptime) 
