library("ballgown")
mCh <- list.files("~/bam_ballgown/stringtie_annotated/", pattern = "mCh.*", full.names = T)
Sxl <- list.files("~/bam_ballgown/stringtie_annotated/", pattern = "Sxl.*", full.names = T)
BAMF <- list.files("~/bam_ballgown/stringtie_annotated/", pattern = "bamF.*", full.names = T)
BAMM <- list.files("~/bam_ballgown/stringtie_annotated/", pattern = "bamM.*", full.names = T)
BG=ballgown(samples = c(BAMF, BAMM, mCh, Sxl), meas = "all")
RNAi_BG=ballgown(samples = c(mCh, Sxl), meas = "all")
BAM_BG=ballgown(samples = c(BAMF, BAMM), meas = "all")
BAMF_Sxl_BG=ballgown(samples = c(BAMF, Sxl), meas = "all")
BAMM_Sxl_BG=ballgown(samples = c(BAMM, Sxl), meas = "all")
BAMF_mCh_BG=ballgown(samples = c(BAMF, mCh), meas = "all")
BAMM_mCh_BG=ballgown(samples = c(BAMM, mCh), meas = "all")
pData(BG) = data.frame(id=sampleNames(BG), group=rep(c("BAMF","BAMM","mCh","Sxl"), each=3))
pData(BAM_BG) = data.frame(id=sampleNames(BAM_BG), group=rep(c("BAMF","BAMM"), each=3))
pData(RNAi_BG) = data.frame(id=sampleNames(RNAi_BG), group=rep(c("mCh","Sxl"), each=3))
pData(BAMF_Sxl_BG) = data.frame(id=sampleNames(BAMF_Sxl_BG), group=rep(c("BAMF","Sxl"), each=3))
pData(BAMM_Sxl_BG) = data.frame(id=sampleNames(BAMM_Sxl_BG), group=rep(c("BAMM","Sxl"), each=3))
pData(BAMF_mCh_BG) = data.frame(id=sampleNames(BAMF_mCh_BG), group=rep(c("BAMF","mCh"), each=3))
pData(BAMM_mCh_BG) = data.frame(id=sampleNames(BAMM_mCh_BG), group=rep(c("BAMM","mCh"), each=3))

eTOt <- data.frame(indexes(BG)$e2t)
tTOg <- data.frame(indexes(BG)$t2g)
eTOg <- merge(eTOt, tTOg, by = "t_id")
eTOg <- eTOg[,-c(1)]
eTOg <- unique(eTOg)
write.table(tTOg, "~/bam_sxl_analysis/data_tables/t2g_annotated.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(eTOg, "~/bam_sxl_analysis/data_tables/e2g_annotated.txt", quote = F, sep = "\t", row.names = F, col.names = T)

merge_gene_stattest <- function(BG_OBJECT, NAME, gene_pval){
  temp_df <- stattest(BG_OBJECT, feature = "gene", meas='FPKM', covariate = 'group', getFC = T)
  temp_df <- temp_df[,c(2:5)]
  colnames(temp_df)[3] <- paste("pval_", NAME, sep = "")
  colnames(temp_df)[4] <- paste("qval_", NAME, sep = "")
  colnames(temp_df)[2] <- paste("fc_", NAME, sep = "")
  temp_df[temp_df[, 2] == 1, 2] <- 0
  gene_pval <- merge(gene_pval, temp_df, by = "id")
  return(gene_pval)
}
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
BG_gene[is.na(BG_gene)] <- 1

merge_trans_stattest <- function(BG_OBJECT, NAME, trans_pval){
  temp_df <- stattest(BG_OBJECT, feature = "transcript", meas = "FPKM", covariate = "group", getFC = T)
  temp_df <- temp_df[,c(2:5)]
  colnames(temp_df)[3] <- paste("pval_", NAME, sep = "")
  colnames(temp_df)[4] <- paste("qval_", NAME, sep = "")
  colnames(temp_df)[2] <- paste("fc_", NAME, sep = "")
  temp_df[temp_df[, 2] == 1, 2] <- 0
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
BG_trans[is.na(BG_trans)] <- 1

merge_exon_stattest <- function(BG_OBJECT, NAME, exon_pval){
  temp_df <- stattest(BG_OBJECT, feature = "exon", covariate = "group", getFC = T)
  temp_df <- temp_df[,c(2:5)]
  colnames(temp_df)[3] <- paste("pval_", NAME, sep = "")
  colnames(temp_df)[4] <- paste("qval_", NAME, sep = "")
  colnames(temp_df)[2] <- paste("fc_", NAME, sep = "")
  temp_df[temp_df[, 2] == 1, 2] <- 0
  exon_pval <- merge(exon_pval, temp_df, by = "id")
  return(exon_pval)
}
BG_exon = stattest(BG, feature = "exon", covariate = "group")
BG_exon <- BG_exon[,c(2,3,4)]
colnames(BG_exon)[2] <- "pval_all"
colnames(BG_exon)[3] <- "qval_all"
BG_exon <- merge_exon_stattest(BAM_BG, "BAM", BG_exon)
BG_exon <- merge_exon_stattest(RNAi_BG, "RNAi", BG_exon)
BG_exon <- merge_exon_stattest(BAMF_Sxl_BG, "BAMF_Sxl", BG_exon)
BG_exon <- merge_exon_stattest(BAMM_Sxl_BG, "BAMM_Sxl", BG_exon)
BG_exon <- merge_exon_stattest(BAMF_mCh_BG, "BAMF_mCh", BG_exon)
BG_exon <- merge_exon_stattest(BAMM_mCh_BG, "BAMM_mCh", BG_exon)
BG_exon[is.na(BG_exon)] <- 1

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

exon_exp_BG <- data.frame(eexpr(BG, 'all'))
colnames(exon_exp_BG)[1] <- "id"
exon_exp_BG$bamF_count_mean <- rowMeans(exon_exp_BG[,c(6, 13, 20)])
exon_exp_BG$bamM_count_mean <- rowMeans(exon_exp_BG[,c(27, 34, 41)])
exon_exp_BG$mCh_count_mean <- rowMeans(exon_exp_BG[,c(48, 55, 62)])
exon_exp_BG$Sxl_count_mean <- rowMeans(exon_exp_BG[,c(69, 76, 83)])
exon_exp_table <- exon_exp_BG
exon_exp_table <- merge(exon_exp_table, BG_exon, by = "id")

write.table(gene_exp_table, "~/bam_sxl_analysis/data_tables/raw_gene_table_annotated.txt", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(trans_exp_table, "~/bam_sxl_analysis/data_tables/raw_trans_table_annotated.txt", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(exon_exp_table, "~/bam_sxl_analysis/data_tables/raw_exon_table_annotated.txt", row.names = F, col.names = T, sep = "\t", quote = F)


