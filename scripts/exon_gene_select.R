exon <- read.delim("~/bam_sxl_analysis/data_tables/raw_exon_table_annotated.txt", header = T)
exon <- exon[,c(1:5, 90:93, 96, 97, 99,100)]
gene <- read.delim("~/bam_sxl_analysis/data_tables/raw_gene_table_annotated.txt", header = T)
gene <- gene[,c(1,14:17, 20, 21, 23, 24)]



pval <- .5
pos_fc <- 1
neg_fc <- 1
count <- 10
fpkm <- .5

exon_pos_fc <- exon[exon$fc_RNAi > pos_fc & exon$pval_RNAi < pval,]
exon_pos_fc <- exon_pos_fc[exon_pos_fc$Sxl_count_mean > count | exon_pos_fc$mCh_count_mean > count,]
exon_neg_fc <- exon[exon$fc_RNAi < neg_fc & exon$pval_RNAi < pval,]
exon_neg_fc <- exon_neg_fc[exon_neg_fc$Sxl_count_mean > count | exon_neg_fc$mCh_count_mean > count,]

gene_pos_fc <- gene[gene$fc_RNAi > pos_fc & gene$pval_RNAi < pval,]
gene_pos_fc <- gene_pos_fc[gene_pos_fc$Sxl_FPKM > fpkm | gene_pos_fc$mCh_FPKM > fpkm, ]
gene_neg_fc <- gene[gene$fc_RNAi < neg_fc & gene$pval_RNAi < pval,]
gene_neg_fc <- gene_neg_fc[gene_neg_fc$Sxl_FPKM > fpkm | gene_neg_fc$mCh_FPKM > fpkm, ]

write.table(exon_pos_fc, "~/bam_sxl_analysis/data_tables/exon_pos_fc_annot.txt", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(exon_neg_fc, "~/bam_sxl_analysis/data_tables/exon_neg_fc_annot.txt", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(gene_pos_fc, "~/bam_sxl_analysis/data_tables/gene_pos_fc_annot.txt", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(gene_neg_fc, "~/bam_sxl_analysis/data_tables/gene_neg_fc_annot.txt", row.names = F, col.names = T, quote = F, sep = "\t")









