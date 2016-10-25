library("ggplot2")
wor_dir <- "/Users/cmdb/sxl_data/graphs/"
gene <- read.delim("~/sxl_data/data_tables/raw_gene_table_annotated.txt", header = T)
gene <- gene[-c(1476, 1468, 286, 7184, 8765),] #FBgn0013686, FBgn0013676, FBgn0001225, FBgn0036091, FBgn0038395
gene_cor <- gene[,c(2:13)]
gene_plot <- log2(gene[,c(2:13)] + 1)
for (i in 1:12){
  for (ii in 2:12) {
    if (ii <= i){next}
    temp_cor <- round(cor(gene_cor[,i], gene_cor[,ii]), 3)
    temp_name <- paste(paste(colnames(gene_plot)[i], colnames(gene_plot)[ii], sep = "_"), ".png", sep = "")
    png(paste(wor_dir, temp_name, sep = ""), width = 6, height = 3.5, units = "in", res = 1200, pointsize = 5)
    print(ggplot(gene_plot, aes(x = gene_plot[,i], y = gene_plot[,ii])) + geom_point(alpha = .15) + geom_smooth(method = "lm") + ggtitle("Sample FPKM Gene Correlation") + labs(x=colnames(gene_plot)[i], y = colnames(gene_plot)[ii]) + annotate("text", x=2,y=8.5, label = paste("r=", temp_cor, sep = "")))
    dev.off()
  }
}

png(paste(wor_dir, "bamM_FPKM-bamF_FPKM.png", sep = ""), width = 6, height = 3.5, units = "in", res = 1200, pointsize = 5)
ggplot(gene, aes(y=log2((bamM_FPKM + 1) / (bamF_FPKM + 1)), x = log2((bamM_FPKM + 1) * (bamF_FPKM + 1)))) + geom_point(aes(color = pval_BAM), alpha = .25) + geom_smooth(color = "red", span = .2) + scale_color_gradient2(name="p_value", low = "red", mid = "gray", high = "blue", midpoint = .5) + ggtitle("Gene FPKM") + labs(x="log2 mean expression", y = "<-BamF log2foldchange BamM->")
dev.off()
png(paste(wor_dir, "bamM_FPKM-mCh_FPKM.png", sep = ""), width = 6, height = 3.5, units = "in", res = 1200, pointsize = 5)
ggplot(gene, aes(y=log2((bamM_FPKM + 1) / (mCh_FPKM + 1)), x = log2((mCh_FPKM + 1) * (bamM_FPKM + 1)))) + geom_point(aes(color = pval_BAMM_mCh), alpha = .25) + geom_smooth(color = "red", span = .2) + scale_color_gradient2(name="p_value", low = "red", mid = "gray", high = "blue", midpoint = .5) + ggtitle("Gene FPKM") + labs(x="log2 mean expression", y = "<-mCh log2foldchange BamM->")
dev.off()
png(paste(wor_dir, "bamF_FPKM-mCh_FPKM.png", sep = ""), width = 6, height = 3.5, units = "in", res = 1200, pointsize = 5)
ggplot(gene, aes(y=log2((bamF_FPKM + 1) / (mCh_FPKM + 1)), x = log2((mCh_FPKM + 1) * (bamF_FPKM + 1)))) + geom_point(aes(color = pval_BAMF_mCh), alpha = .25) + geom_smooth(color = "red", span = .2) + scale_color_gradient2(name="p_value", low = "red", mid = "gray", high = "blue", midpoint = .5) + ggtitle("Gene FPKM") + labs(x="log2 mean expression", y = "<-mCh log2foldchange BamF->")
dev.off()
png(paste(wor_dir, "bamM_FPKM-Sxl_FPKM.png", sep = ""), width = 6, height = 3.5, units = "in", res = 1200, pointsize = 5)
ggplot(gene, aes(y=log2((bamM_FPKM + 1) / (Sxl_FPKM + 1)), x = log2((Sxl_FPKM + 1) * (bamM_FPKM + 1)))) + geom_point(aes(color = pval_BAMM_Sxl), alpha = .25) + geom_smooth(color = "red", span = .2) + scale_color_gradient2(name="p_value", low = "red", mid = "gray", high = "blue", midpoint = .5) + ggtitle("Gene FPKM") + labs(x="log2 mean expression", y = "<-Sxl log2foldchange BamM->")
dev.off()
png(paste(wor_dir, "bamF_FPKM-Sxl_FPKM.png", sep = ""), width = 6, height = 3.5, units = "in", res = 1200, pointsize = 5)
ggplot(gene, aes(y=log2((bamF_FPKM + 1) / (Sxl_FPKM + 1)), x = log2((Sxl_FPKM + 1) * (bamF_FPKM + 1)))) + geom_point(aes(color = pval_BAMF_Sxl), alpha = .25) + geom_smooth(color = "red", span = .2) + scale_color_gradient2(name="p_value", low = "red", mid = "gray", high = "blue", midpoint = .5) + ggtitle("Gene FPKM") + labs(x="log2 mean expression", y = "<-Sxl log2foldchange BamF->")
dev.off()
png(paste(wor_dir, "mCh_FPKM-Sxl_FPKM.png", sep = ""), width = 6, height = 3.5, units = "in", res = 1200, pointsize = 5)
ggplot(gene, aes(y=log2((mCh_FPKM + 1) / (Sxl_FPKM + 1)), x = log2((Sxl_FPKM + 1) * (mCh_FPKM + 1)))) + geom_point(aes(color = pval_RNAi), alpha = .25) + geom_smooth(color = "red", span = .2) + scale_color_gradient2(name="p_value", low = "red", mid = "gray", high = "blue", midpoint = .5) + ggtitle("Gene FPKM") + labs(x="log2 mean expression", y = "<-Sxl log2foldchange mCh->")
dev.off()

a <- gene[gene$id == "FBgn0264270",]
bamf <- data.frame(rowMeans(a[,2:4]))
bamfsd <- sd(a[,2:4])
bamf <- bamf[1,1]
bamm <- data.frame(rowMeans(a[,5:7]))
bammsd <- sd(a[,5:7])
bamm <- bamm[1,1]
mch <- data.frame(rowMeans(a[,8:10]))
mchsd <- sd(a[,8:10])
mch <- mch[1,1]
sxl <- data.frame(rowMeans(a[,11:13]))
sxlsd <- sd(a[,11:13])
sxl <- sxl[1,1]
temp_means <- c(bamf, bamm, mch, sxl)
temp_names <- c("bamf", "bamm", "mch", "sxl")
df <- data.frame(temp_means)
rownames(df) <- temp_names
colnames(df)[1] <- "means"
df$sd <- c(bamfsd, bammsd, mchsd, sxlsd)
df$ymax <- df$means + df$sd
df$ymin <- df$means - df$sd
df$names <- rownames(df)
limits <- aes(ymax = df$ymax, ymin = df$ymin)
png("/Users/cmdb/sxl_data/graphs/sxl_expression.png", width = 6, height = 3.5, units = "in", res = 1200, pointsize = 5)
ggplot(df, aes(x = names, y = means)) + geom_bar( stat = "identity") + geom_errorbar(limits, width = .15) + ggtitle("Sxl Expression") + labs(x = "Genotypes", y = "FPKM Mean")
dev.off()



a <- gene[gene$mCh_FPKM > 0 & gene$Sxl_FPKM > 0,]
hist(log2(a$Sxl_FPKM + 1), breaks = 100, ylim = c(0,2000), xlim = c(0,15))
hist(log2(a$mCh_FPKM + 1), breaks = 100, ylim = c(0,2000), xlim = c(0,15))
a <- gene[gene$bamF_FPKM > 0 & gene$bamM_FPKM > 0,]
hist(log2(a$bamF_FPKM + 1), breaks = 100, ylim = c(0,2000), xlim = c(0,15))
hist(log2(a$bamM_FPKM + 1), breaks = 200, ylim = c(0,2000), xlim = c(0,15))


df <- gene[,c(14:17)]
rownames(df) <- gene$id
df <- log2(df[rowSums(df[,1:4]) > 0,] + 1)
#df <- df[,-c(5)]
#mydata <- scale(df)
mydata <- df
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata, centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
fit <- kmeans(mydata, 4)
aggregate(mydata,by=list(fit$cluster), FUN=mean)
kmeans_data <- data.frame(mydata, fit$cluster)


df <- gene[,c(14,16,17,15)]
rownames(df) <- gene$id
df <- log2(df[rowSums(df[,1:4]) > 4,] + 1)
#df <- df[,-c(5)]
#df <- df[df$bamF_FPKM < 100 | df$bamM_FPKM < 100 | df$Sxl_FPKM < 100 | df$mCh_FPKM < 100,]
#mydata <- scale(df)
mydata <- df
library("reshape")
library("scales")
ord <- hclust( dist(mydata, method = "euclidean"), method = "ward.D" )$order
pd <- as.data.frame(mydata)
pd$gene <- rownames(pd)
pd.m <- melt(pd, id.vars = "gene", variable = "genotype" )
pd.m$genotype <- factor( pd.m$genotype, levels = colnames(mydata), labels = colnames(mydata) )
pd.m$gene <- factor( pd.m$gene, levels = rownames(mydata)[ord],  labels = seq_along(rownames(mydata)) )
temp <- data.frame(quantile(pd.m$value))
low <- round(temp[1,1], 3)
mid <- round(temp[3,1] + .001, 3)
high <- round(temp[5,1] + .001, 3)
png("/Users/cmdb/sxl_data/graphs/clust_gene_ward.d.png", width = 6, height = 3.5, units = "in", res = 1200, pointsize = 5)
ggplot(pd.m, aes(gene, genotype)) + geom_tile(aes(fill = value)) + scale_fill_gradientn(colors = c("blue", "grey", "red"), values = rescale(c(low, mid, high)))
dev.off()

ggplot(pd.m, aes(gene, genotype)) + geom_tile(aes(fill = value)) + scale_fill_gradient2(low=muted("blue"), midpoint = 6.25, mid = "grey", high = muted("red"))
scale_fill_gradient2(colors = c("blue", "red", "red"), values = rescale(c(-.2, .3, 50)))

################################################################################
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
exon <- read.delim("~/sxl_data/data_tables/raw_exon_table_annotated.txt", header = T)
sxl <- exon[exon$id == "73916",]
sxl_read <- sxl[c(6,13,20,27,34,41,48,55,62,69,76,83)]
sxl_cov <- sxl[c(9,16,23,30,37,44,51,58,65,72,79,86)]
a <- sxl_read

bamf <- data.frame(rowMeans(a[,1:3]))
bamfsd <- sd(a[,1:3])
bamf <- bamf[1,1]
bamm <- data.frame(rowMeans(a[,4:6]))
bammsd <- sd(a[,4:6])
bamm <- bamm[1,1]
mch <- data.frame(rowMeans(a[,7:9]))
mchsd <- sd(a[,7:9])
mch <- mch[1,1]
sxl <- data.frame(rowMeans(a[,10:12]))
sxlsd <- sd(a[,10:12])
sxl <- sxl[1,1]
temp_means <- c(bamf, bamm, mch, sxl)
temp_names <- c("bamf", "bamm", "mch", "sxl")
df <- data.frame(temp_means)
rownames(df) <- temp_names
colnames(df)[1] <- "means"
df$sd <- c(bamfsd, bammsd, mchsd, sxlsd)
df$ymax <- df$means + df$sd
df$ymin <- df$means - df$sd
df$names <- rownames(df)
limits <- aes(ymax = df$ymax, ymin = df$ymin)
png("/Users/cmdb/sxl_data/graphs/sxl_male_exon_count.png", width = 6, height = 3.5, units = "in", res = 1200, pointsize = 5)
ggplot(df, aes(x = names, y = means)) + geom_bar( stat = "identity") + geom_errorbar(limits, width = .15) + ggtitle("Sxl Male Exon Count") + labs(x = "Genotypes", y = "Exon Count")
dev.off()




################################################################################################ 

################################################################################################ 

