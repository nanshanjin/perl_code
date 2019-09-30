require(ggplot2)
require(RColorBrewer)
args <- commandArgs(TRUE)
datafile <- args[1]
data <- read.table(datafile, header=TRUE, sep='\t')
#根据biological_process筛选
df_BP <- subset(data, GO_namespace == 'biological_process')
#order利用索引排序，decreasing为TRUE，表示降序
order_df_BP <- df_BP[order(df_BP$gene_num, decreasing = TRUE),]

order_df_BP$Term <- as.vector(order_df_BP$Term)
order_df_BP$gene_num <- as.numeric(order_df_BP$gene_num)
myLabel <- paste(order_df_BP$Term, " (", round(order_df_BP$gene_num / sum(order_df_BP$gene_num) * 100, 2), "%)", sep = "")
myLabel <- myLabel[order(myLabel)]
newpalette <- colorRampPalette(brewer.pal(12,"Set3"))(length(order_df_BP$Term))
p<-ggplot(order_df_BP, aes(x = '', y = gene_num, fill = Term)) + geom_bar(stat = 'identity', width = 1) + coord_polar(theta = 'y') + labs(x = "", y = "", title = "Biological Process Pie Chart") + theme_bw() + theme(axis.ticks = element_blank(), legend.title = element_blank(), legend.position = "right", axis.text.x = element_blank(), panel.border=element_blank(), plot.title = element_text(hjust = 0.5), panel.grid=element_blank()) + scale_fill_discrete(breaks = order_df_BP$gene_num, labels = myLabel) + scale_fill_manual(values = newpalette,labels = myLabel)
ggsave("BP_level.png", plot = p, width = 20, height = 10)
ggsave("BP_level.pdf", plot = p, width = 20, height = 10)
df_CC <- subset(data, GO_namespace == 'cellular_component')
order_df_CC <- df_CC[order(df_CC$gene_num, decreasing = TRUE),]
order_df_CC$Term <- as.vector(order_df_CC$Term)
order_df_CC$gene_num <- as.numeric(order_df_CC$gene_num)
myLabel <- paste(order_df_CC$Term, " (", round(order_df_CC$gene_num / sum(order_df_CC$gene_num) * 100, 2), "%)", sep = "")
myLabel <- myLabel[order(myLabel)]
newpalette <- colorRampPalette(brewer.pal(12,"Set3"))(length(order_df_CC$Term))
p<-ggplot(order_df_CC, aes(x = '', y = gene_num, fill = Term)) + geom_bar(stat = 'identity', width = 1) + coord_polar(theta = 'y') + labs(x = "", y = "", title = "Cellular Component Pie Chart") + theme_bw() + theme(axis.ticks = element_blank(), legend.title = element_blank(), legend.position = "right", axis.text.x = element_blank(), panel.border=element_blank(), plot.title = element_text(hjust = 0.5), panel.grid=element_blank()) + scale_fill_discrete(breaks = order_df_CC$gene_num, labels = myLabel) + scale_fill_manual(values = newpalette,labels = myLabel)
ggsave("CC_level.png", plot = p, width = 20, height = 10)
ggsave("CC_level.pdf", plot = p, width = 20, height = 10)
df_MF <- subset(data, GO_namespace == 'molecular_function')
order_df_MF <- df_MF[order(df_MF$gene_num, decreasing = TRUE),]
order_df_MF$Term <- as.vector(order_df_MF$Term)
order_df_MF$gene_num <- as.numeric(order_df_MF$gene_num)
myLabel <- paste(order_df_MF$Term, " (", round(order_df_MF$gene_num / sum(order_df_MF$gene_num) * 100, 2), "%)", sep = "")
myLabel <- myLabel[order(myLabel)]
newpalette <- colorRampPalette(brewer.pal(12,"Set3"))(length(order_df_MF$Term))
p<-ggplot(order_df_MF, aes(x = '', y = gene_num, fill = Term)) + geom_bar(stat = 'identity', width = 1) + coord_polar(theta = 'y') + labs(x = "", y = "", title = "Molecular Function Pie Chart") + theme_bw() + theme(axis.ticks = element_blank(), legend.title = element_blank(), legend.position = "right", axis.text.x = element_blank(), panel.border=element_blank(), plot.title = element_text(hjust = 0.5), panel.grid=element_blank()) + scale_fill_discrete(breaks = order_df_MF$gene_num, labels = myLabel) + scale_fill_manual(values = newpalette,labels = myLabel)
ggsave("MF_level.png", plot = p, width = 20, height = 10)
ggsave("MF_level.pdf", plot = p, width = 20, height = 10)
