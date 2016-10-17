library("ballgown")
library("ggplot2")
RNAi_BG=ballgown(samples = list.files("~/bam_ballgown", pattern = ".*RNAi.*", full.names = T), meas = "all")
pData(RNAi_BG) = data.frame(id=sampleNames(RNAi_BG), group=rep(c("mCh","Sxl"), each=3))


eTOt <- data.frame(indexes(RNAi_BG)$e2t)
tTOg <- data.frame(indexes(RNAi_BG)$t2g)
eTOg <- merge(eTOt, tTOg, by = "t_id")
eTOg <- eTOg[,-c(1)]
eTOg <- unique(eTOg)


e <- data.frame(eexpr(RNAi_BG, 'all'))
e$mCh_count_mean <- rowMeans(e[,c(6, 13, 20)])
e$Sxl_count_mean <- rowMeans(e[,c(27, 34, 41)])
e_stats <- stattest(RNAi_BG, feature = 'exon', covariate = 'group', meas = 'rcount', getFC = T)
e_stats <- e_stats[,-c(1)]
colnames(e_stats)[1] <- "e_id"
e_data <- merge(e, e_stats, by = "e_id")
e_data[is.na(e_data)] <- 1
e_data$fc[e_data$fc == 1] <- 0
write.table(e_sig, "~/bam_sxl_analysis/data_tables/exons.txt", row.names = F, col.names = T, sep = "\t", quote = F)
e_sig <- e_data[e_data$pval < 0.05 & e_data$fc > 1, ]
e_sig <- e_sig[e_sig$Sxl_count_mean > 10 | e_sig$mCh_count_mean > 10,]
write.table(e_sig, "~/bam_sxl_analysis/data_tables/sig_exons.txt", row.names = F, col.names = T, sep = "\t", quote = F)



################################################################################################

targets <- read.table("~/bam_sxl_analysis/data_tables/sxl_targets.txt", sep = "\t", header = F)
colnames(targets) <- c("e_id", "location", "position", "sequence")
targets <- merge(targets, eTOg, by = "e_id")
targets <- merge(targets, e_sig, by = "e_id")
write.table(targets, "~/bam_sxl_analysis/data_tables/sxl_targets_annotated.txt", row.names = F, col.names = T, sep = "\t", quote = F)
################################################################################################


################################################################################################


BG_gene = stattest(BG, feature = "gene", meas='FPKM', covariate = 'group')
BG_gene <- BG_gene[,c(2,3,4)]
colnames(BG_gene)[2] <- "pval_all"
colnames(BG_gene)[3] <- "qval_all"
BG_gene <- merge_gene_stattest(BAM_BG, "BAM", BG_gene)
BG_gene <- merge_gene_stattest(RNAi_BG, "RNAi", BG_gene)
BG_gene <- merge_gene_stattest(BAMF_Sxl_BG, "BAMF_Sxl", BG_gene)
BG_gene <- merge_gene_stattest(BAMM_Sxl_BG, "BAMM_Sxl", BG_gene)
BG_gene <- merge_gene_stattest(BAMF_mCh_BG, "BAMF_mCh", BG_gene)
BG_gene <- merge_gene_stattest(BAMM_mCh_BG, "BAMM_mCh", BG_gene)

merge_trans_stattest <- function(BG_OBJECT, NAME, trans_pval){
  temp_df <- stattest(BG_OBJECT, feature = "transcript", meas = "FPKM", covariate = "group")
  temp_df <- temp_df[,c(2,3,4)]
  colnames(temp_df)[2] <- paste("pval_", NAME, sep = "")
  colnames(temp_df)[3] <- paste("qval_", NAME, sep = "")
  trans_pval <- merge(trans_pval, temp_df, by = "id")
  return(trans_pval)
}
BG_trans = stattest(BG, feature = "transcript", meas = "FPKM", covariate = "group")
BG_trans <- BG_trans[,c(2,3,4)]
colnames(BG_trans)[2] <- "pval_all"
colnames(BG_trans)[3] <- "qval_all"
BG_trans <- merge_trans_stattest(BAM_BG, "BAM", BG_trans)
BG_trans <- merge_trans_stattest(RNAi_BG, "RNAi", BG_trans)
BG_trans <- merge_trans_stattest(BAMF_Sxl_BG, "BAMF_Sxl", BG_trans)
BG_trans <- merge_trans_stattest(BAMM_Sxl_BG, "BAMM_Sxl", BG_trans)
BG_trans <- merge_trans_stattest(BAMF_mCh_BG, "BAMF_mCh", BG_trans)
BG_trans <- merge_trans_stattest(BAMM_mCh_BG, "BAMM_mCh", BG_trans)
BG_gene[is.na(BG_gene)] <- 1
BG_trans[is.na(BG_trans)] <- 1

gene_exp_BG = data.frame(gexpr(BG))
gene_exp_BG$bamF_FPKM <- rowMeans(gene_exp_BG[, c(1:3)])
gene_exp_BG$bamM_FPKM <- rowMeans(gene_exp_BG[, c(4:6)])
gene_exp_BG$mCh_FPKM <- rowMeans(gene_exp_BG[, c(7:9)])
gene_exp_BG$Sxl_FPKM <- rowMeans(gene_exp_BG[, c(10:12)])
gene_exp_BG$id <- rownames(gene_exp_BG)
rownames(gene_exp_BG) <- NULL
gene_exp_table <- gene_exp_BG[,c(17, 1:16)]
gene_exp_table <- merge(gene_exp_table, BG_gene, by = "id")

trans_exp_BG = data.frame(texpr(BG, meas = "all"))
colnames(trans_exp_BG)[1] <- "id"
trans_exp_BG$bamF_FPKM <- rowMeans(trans_exp_BG[, c(12,14,16)])
trans_exp_BG$bamM_FPKM <- rowMeans(trans_exp_BG[, c(18,20,22)])
trans_exp_BG$mCh_FPKM <- rowMeans(trans_exp_BG[, c(24,26,28)])
trans_exp_BG$Sxl_FPKM <- rowMeans(trans_exp_BG[, c(30,32,34)])
trans_exp_table <- trans_exp_BG
trans_exp_table <- merge(trans_exp_table, BG_trans, by = "id")

BG_exon <- stattest(BG, feature = "exon", covariate = "group")


gene_exp_table_modified <- gene_exp_table[-c(1476, 7184),] #FBgn0013686, FBgn0036091
trans_exp_table_modified <- trans_exp_table[-c(35069, 16877),] #FBgn0013686, FBgn0036091
################################################################################################ 
write.table(gene_exp_table_modified, "~/bam_sxl_analysis/data_tables/raw_gene_table.txt", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(trans_exp_table_modified, "~/bam_sxl_analysis/data_tables/raw_trans_table.txt", row.names = F, col.names = T, sep = "\t", quote = F)
################################################################################################ 
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

