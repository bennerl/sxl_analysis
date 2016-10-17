exon <- read.delim("~/sxl_data/data_tables/raw_exon_table_unannotated.txt", header = T)
exon <- exon[,c(1:5, 90:93, 96, 97, 99,100)]

pval <- .5
pos_fc <- 1
neg_fc <- 1
count <- 10

exon_pos_fc <- exon[exon$fc_RNAi > pos_fc & exon$pval_RNAi < pval,]
exon_pos_fc <- exon_pos_fc[exon_pos_fc$Sxl_count_mean > count | exon_pos_fc$mCh_count_mean > count,]
exon_neg_fc <- exon[exon$fc_RNAi < neg_fc & exon$pval_RNAi < pval,]
exon_neg_fc <- exon_neg_fc[exon_neg_fc$Sxl_count_mean > count | exon_neg_fc$mCh_count_mean > count,]

write.table(exon_pos_fc, "~/sxl_data/data_tables/exon_pos_fc_unannot.txt", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(exon_neg_fc, "~/sxl_data/data_tables/exon_neg_fc_unannot.txt", row.names = F, col.names = T, quote = F, sep = "\t")









