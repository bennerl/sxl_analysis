names <- read.delim("~/annotations/Dmel_gene_names.txt", header = F)
names <- names[,c(1,2)]
colnames(names) <- c("name","id")
id <- read.delim("~/sxl_data/data_tables/e2g_annotated.txt", header = T)
colnames(id)[2] <- "id"
tran_id <- read.delim("~/sxl_data/data_tables/t2g_annotated.txt", header = T)
colnames(tran_id)[2] <- "id"
exon_pos <- read.delim("~/sxl_data/data_tables/exon_pos_fc_annot.txt", header = T)
colnames(exon_pos)[1] <- "e_id"
exon_pos <- merge(id, exon_pos, by = "e_id")
exon_pos <- merge(names, exon_pos, by = "id")
write.table(exon_pos, "~/sxl_data/data_tables/name_exon_pos_fc_annot.txt", row.names = F, col.names = T, sep = "\t", quote = F)
exon_neg <- read.delim("~/sxl_data/data_tables/exon_neg_fc_annot.txt", header = T)
colnames(exon_neg)[1] <- "e_id"
exon_neg <- merge(id, exon_neg, by = "e_id")
exon_neg <- merge(names, exon_neg, by = "id")
write.table(exon_neg, "~/sxl_data/data_tables/name_exon_neg_fc_annot.txt", row.names = F, col.names = T, sep = "\t", quote = F)
gene_pos <- read.delim("~/sxl_data/data_tables/gene_pos_fc_annot.txt", header = T)
gene_pos <- merge(names, gene_pos, by = "id")
write.table(gene_pos, "~/sxl_data/data_tables/name_gene_pos_fc_annot.txt", row.names = F, col.names = T, sep = "\t", quote = F)
gene_neg <- read.delim("~/sxl_data/data_tables/gene_neg_fc_annot.txt", header = T)
gene_neg <- merge(names, gene_neg, by = "id")
write.table(gene_neg, "~/sxl_data/data_tables/name_gene_neg_fc_annot.txt", row.names = F, col.names = T, sep = "\t", quote = F)
raw_exon <- read.delim("~/sxl_data/data_tables/raw_exon_table_annotated.txt", header = T)
colnames(raw_exon)[1] <- "e_id"
raw_exon <- merge(id, raw_exon, by = "e_id")
raw_exon <- merge(names, raw_exon, by = "id")
write.table(raw_exon, "~/sxl_data/data_tables/name_raw_exon_table_annotated.txt", row.names = F, col.names = T, sep = "\t", quote = F)
raw_gene <- read.delim("~/sxl_data/data_tables/raw_gene_table_annotated.txt", header = T)
raw_gene <- merge(names, raw_gene, by = "id")
write.table(raw_gene, "~/sxl_data/data_tables/name_raw_gene_table_annotated.txt", row.names = F, col.names = T, sep = "\t", quote = F)
raw_trans <- read.delim("~/sxl_data/data_tables/raw_trans_table_annotated.txt", header = T)
colnames(raw_trans)[1] <- "t_id"
raw_trans <- merge(tran_id, raw_trans, by = "t_id")
raw_trans <- merge(names, raw_trans, by = "id")
write.table(raw_trans, "~/sxl_data/data_tables/name_raw_trans_table_annotated.txt", row.names = F, col.names = T, sep = "\t", quote = F)


