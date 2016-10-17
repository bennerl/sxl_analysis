gene <- gene_exp_table_modified
library("GGally")
A <- log2(gene[,2:17]+1)
plt <- ggpairs(A[,13:16], columnLabels = c("bamF","bamM","mCh","Sxl"), axisLabels = "none", title = "Genotype FPKM")
cor_vals <- ggpairs(gene[,14:17] )
for (i in 1:4){
  for (ii in 2:4){
    if (ii <= i) {next}
    temp_val <- getPlot(cor_vals, i, ii)
    plt <- putPlot(plt, temp_val, i,ii)
  }
}
plt

plt <- ggpairs(A[,1:12], columnLabels = c("F1","F2","F3","M1","M2","M3","C1","C2","C3","S1","S2","S3"), axisLabels = "none", title = "Genotype FPKM")
cor_vals <- ggpairs(gene[,2:13])
for (i in 1:12){
  for (ii in 2:12){
    if (ii <= i) {next}
    temp_val <- getPlot(cor_vals, i, ii)
    plt <- putPlot(plt, temp_val, i,ii)
  }
}
plt
################################################################################################ 

ggplot(gene, aes(y=log2((bamM_FPKM + 1) / (bamF_FPKM + 1)), x = log2((bamM_FPKM + 1) * (bamF_FPKM + 1)))) + geom_point(aes(color = pval_BAM), alpha = .25) + geom_smooth(color = "red", span = .2) + scale_color_gradient2(name="p_value", low = "red", mid = "gray", high = "blue", midpoint = .5) + ggtitle("Gene FPKM") + labs(x="log2 mean expression", y = "<-BamF log2foldchange BamM->")
ggplot(gene, aes(y=log2((bamM_FPKM + 1) / (mCh_FPKM + 1)), x = log2((mCh_FPKM + 1) * (bamM_FPKM + 1)))) + geom_point(aes(color = pval_BAMM_mCh), alpha = .25) + geom_smooth(color = "red", span = .2) + scale_color_gradient2(name="p_value", low = "red", mid = "gray", high = "blue", midpoint = .5) + ggtitle("Gene FPKM") + labs(x="log2 mean expression", y = "<-mCh log2foldchange BamM->")
ggplot(gene, aes(y=log2((bamF_FPKM + 1) / (mCh_FPKM + 1)), x = log2((mCh_FPKM + 1) * (bamF_FPKM + 1)))) + geom_point(aes(color = pval_BAMF_mCh), alpha = .25) + geom_smooth(color = "red", span = .2) + scale_color_gradient2(name="p_value", low = "red", mid = "gray", high = "blue", midpoint = .5) + ggtitle("Gene FPKM") + labs(x="log2 mean expression", y = "<-mCh log2foldchange BamF->")
ggplot(gene, aes(y=log2((bamM_FPKM + 1) / (Sxl_FPKM + 1)), x = log2((Sxl_FPKM + 1) * (bamM_FPKM + 1)))) + geom_point(aes(color = pval_BAMM_Sxl), alpha = .25) + geom_smooth(color = "red", span = .2) + scale_color_gradient2(name="p_value", low = "red", mid = "gray", high = "blue", midpoint = .5) + ggtitle("Gene FPKM") + labs(x="log2 mean expression", y = "<-Sxl log2foldchange BamM->")
ggplot(gene, aes(y=log2((bamF_FPKM + 1) / (Sxl_FPKM + 1)), x = log2((Sxl_FPKM + 1) * (bamF_FPKM + 1)))) + geom_point(aes(color = pval_BAMF_Sxl), alpha = .25) + geom_smooth(color = "red", span = .2) + scale_color_gradient2(name="p_value", low = "red", mid = "gray", high = "blue", midpoint = .5) + ggtitle("Gene FPKM") + labs(x="log2 mean expression", y = "<-Sxl log2foldchange BamF->")
ggplot(gene, aes(y=log2((mCh_FPKM + 1) / (Sxl_FPKM + 1)), x = log2((Sxl_FPKM + 1) * (mCh_FPKM + 1)))) + geom_point(aes(color = pval_RNAi), alpha = .25) + geom_smooth(color = "red", span = .2) + scale_color_gradient2(name="p_value", low = "red", mid = "gray", high = "blue", midpoint = .5) + ggtitle("Gene FPKM") + labs(x="log2 mean expression", y = "<-Sxl log2foldchange mCh->")

################################################################################################ 

a <- gene[gene$id == "FBgn0264270",]
write.table(a, "~/asdfasda.txt", sep = "\t", col.names = F, row.names = T)
a<-read.delim("~/asdfasda.txt", header = F, sep = "\t")
A<-read.delim("~/hgefs.txt", header = T, sep = "\t")
B <- data.frame(t(A))
B$means <- rowMeans(B[,c(1:3)])
for (i in 1:nrow(B)){
  B[i,5] <- B[i,4] + sd(B[i,1:3])
  B[i,6] <- B[i,4] - sd(B[i,1:3])
}
colnames(B)[5] <- "ymax"
colnames(B)[6] <- "ymin"
B$names <- rownames(B)
limits <- aes(ymax = B$ymax, ymin = B$ymin)
ggplot(B, aes(x = names, y = means)) + geom_bar( stat = "identity") + geom_errorbar(limits, width = .15) + ggtitle("Sxl Expression") + labs(x = "Genotypes", y = "FPKM Mean")

################################################################################################ 

