#change motif in target_finder.py and 

"""
change description and tracknames for bed files with motifs
change motif based on paper
look at the distribuition of FPKM and read counts for genes, come up with better rationale to remove genes that are underexrpressed in samples
"""



MOTIF="Sakashita"
RScript /Users/cmdb/sxl_analysis/scripts/exon_gene_select.R #
RScript /Users/cmdb/sxl_analysis/scripts/exon_unannotated_gene_select.R #

rm /Users/cmdb/sxl_data/data_tables/gene_*_annot.chr.txt
for i in /Users/cmdb/sxl_data/data_tables/gene_*_annot.txt
do
	/Users/cmdb/sxl_analysis/scripts/add_chr_to_BG_gene_list.py $i > temp.txt
	mv temp.txt ${i%txt}chr.txt
done

rm ~/sxl_data/bed/${MOTIF}_exon*annot.bed
~/sxl_analysis/scripts/target_finder.py exon bed ~/sxl_data/data_tables/exon_pos_fc_annot.txt ~/sxl_data/data_tables/e2g_annotated.txt > ~/sxl_data/bed/${MOTIF}_exon_pos_annot.bed
~/sxl_analysis/scripts/target_finder.py exon bed ~/sxl_data/data_tables/exon_neg_fc_annot.txt ~/sxl_data/data_tables/e2g_annotated.txt > ~/sxl_data/bed/${MOTIF}_exon_neg_annot.bed

echo "track name=\"${MOTIF}ExonAnnot\" description=\"${MOTIF} Exon Annotated\" visibility=2 itemRgb=\"On\"" > ~/sxl_data/bed/${MOTIF}_exon_annot.bed
for i in ~/sxl_data/bed/${MOTIF}_exon_*_annot.bed
do
	FILE=${i#/Users/cmdb/sxl_data/bed/${MOTIF}_exon_}
	FILE=${FILE%_annot.bed}
	if [ $FILE == "pos" ]
	then
		cat $i | sed 's|0,0,0|0,178,46|g'| sort | uniq >> ~/sxl_data/bed/${MOTIF}_exon_annot.bed
	else
		cat $i | sed 's|0,0,0|178,0,0|g'| sort | uniq >> ~/sxl_data/bed/${MOTIF}_exon_annot.bed
	fi
done

rm ~/sxl_data/bed/${MOTIF}_gene*annot.bed
~/sxl_analysis/scripts/target_finder.py gene bed ~/sxl_data/data_tables/gene_pos_fc_annot.chr.txt ~/annotations/r6.12_dmel_noERCC.gtf > ~/sxl_data/bed/${MOTIF}_gene_pos_annot.bed
~/sxl_analysis/scripts/target_finder.py gene bed ~/sxl_data/data_tables/gene_neg_fc_annot.chr.txt ~/annotations/r6.12_dmel_noERCC.gtf > ~/sxl_data/bed/${MOTIF}_gene_neg_annot.bed

echo "track name=\"${MOTIF}GeneAnnot\" description=\"${MOTIF} Gene Annotated\" visibility=2 itemRgb=\"On\"" > ~/sxl_data/bed/${MOTIF}_gene_annot.bed
for i in ~/sxl_data/bed/${MOTIF}_gene_*_annot.bed
do
	FILE=${i#/Users/cmdb/sxl_data/bed/${MOTIF}_gene_}
	FILE=${FILE%_annot.bed}
	if [ $FILE == "pos" ]
	then
		cat $i | sed 's|0,0,0|0,178,46|g'| sort | uniq >> ~/sxl_data/bed/${MOTIF}_gene_annot.bed
	else
		cat $i | sed 's|0,0,0|178,0,0|g'| sort | uniq >> ~/sxl_data/bed/${MOTIF}_gene_annot.bed
	fi
done

rm ~/sxl_data/data_tables/novel_exons.txt
~/sxl_analysis/scripts/novel_exon_find.py gtf > ~/sxl_data/data_tables/novel_exons.txt
rm ~/sxl_data/data_tables/novel_exon_*_fc_unannot.txt
~/sxl_analysis/scripts/novel_exon_find.py merge ~/sxl_data/data_tables/exon_pos_fc_unannot.txt > ~/sxl_data/data_tables/novel_exon_pos_fc_unannot.txt
~/sxl_analysis/scripts/novel_exon_find.py merge ~/sxl_data/data_tables/exon_neg_fc_unannot.txt > ~/sxl_data/data_tables/novel_exon_neg_fc_unannot.txt
rm ~/sxl_data/bed/${MOTIF}_novel_exon*unannot.bed
~/sxl_analysis/scripts/target_novel_exon.py bed ~/sxl_data/data_tables/novel_exon_pos_fc_unannot.txt > ~/sxl_data/bed/${MOTIF}_novel_exon_pos_unannot.bed
~/sxl_analysis/scripts/target_novel_exon.py bed ~/sxl_data/data_tables/novel_exon_neg_fc_unannot.txt > ~/sxl_data/bed/${MOTIF}_novel_exon_neg_unannot.bed

echo "track name=\"${MOTIF}ExonUnannot\" description=\"${MOTIF} Exon Unannotated\" visibility=2 itemRgb=\"On\"" > ~/sxl_data/bed/${MOTIF}_novel_exon_unannot.bed
for i in ~/sxl_data/bed/${MOTIF}_novel_exon_*_unannot.bed
do
	FILE=${i#/Users/cmdb/sxl_data/bed/${MOTIF}_novel_exon_}
	FILE=${FILE%_unannot.bed}
	if [ $FILE == "pos" ]
	then
		cat $i | sed 's|0,0,0|0,178,46|g'| sort | uniq >> ~/sxl_data/bed/${MOTIF}_novel_exon_unannot.bed
	else
		cat $i | sed 's|0,0,0|178,0,0|g'| sort | uniq >> ~/sxl_data/bed/${MOTIF}_novel_exon_unannot.bed
	fi
done

rm ~/sxl_data/data_tables/${MOTIF}_exon_*_annot_stat*.txt
~/sxl_analysis/scripts/target_finder.py exon stats ~/sxl_data/data_tables/exon_pos_fc_annot.txt ~/sxl_data/data_tables/e2g_annotated.txt > ~/sxl_data/data_tables/${MOTIF}_exon_pos_annot_stats.txt
~/sxl_analysis/scripts/target_finder.py exon stats ~/sxl_data/data_tables/exon_neg_fc_annot.txt ~/sxl_data/data_tables/e2g_annotated.txt > ~/sxl_data/data_tables/${MOTIF}_exon_neg_annot_stats.txt

for i in ~/sxl_data/data_tables/${MOTIF}_exon_*stat*.txt
do
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "Gene_id" "Gene" "Motif" "FC_RNAi" "pVal_RNAi" "FC_bam" "pVal_bam" "mCh_FPKM" "Sxl_FPKM" "bamF_FPKM" "bamM_FPKM" "Chr" "Start" "End" "Strand" "Exon_start" "Exon_end" "Prime" "e_id" > temp.txt
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "Gene_id" "Gene" "Motif" "FC_RNAi" "pVal_RNAi" "FC_bam" "pVal_bam" "mCh_FPKM" "Sxl_FPKM" "bamF_FPKM" "bamM_FPKM" "Chr" "Start" "End" "Strand" "exon_start" "exon_end" "Prime" "e_id" > temp1.txt
	FC=${i#/Users/cmdb/sxl_data/data_tables/${MOTIF}_exon_}
	FC=${FC%_annot_stats.txt}
	if [ $FC == "pos" ]
	then
		cat $i | sort -nrk4 | uniq >> temp1.txt
	else
		cat $i | sort -nk4 | uniq >> temp1.txt
	fi
	mv temp1.txt ${i%txt}sort_uniq.txt
	cat $i >> temp.txt; rm $i; mv temp.txt $i
done

rm ~/sxl_data/data_tables/${MOTIF}_gene_*_annot_stat*.txt
~/sxl_analysis/scripts/target_finder.py gene stats ~/sxl_data/data_tables/gene_pos_fc_annot.chr.txt ~/annotations/r6.12_dmel_noERCC.gtf > ~/sxl_data/data_tables/${MOTIF}_gene_pos_annot_stats.txt
~/sxl_analysis/scripts/target_finder.py gene stats ~/sxl_data/data_tables/gene_neg_fc_annot.chr.txt ~/annotations/r6.12_dmel_noERCC.gtf > ~/sxl_data/data_tables/${MOTIF}_gene_neg_annot_stats.txt

for i in ~/sxl_data/data_tables/${MOTIF}_gene*stat*.txt
do
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "Gene_id" "Gene" "Motif" "FC_RNAi" "pVal_RNAi" "FC_bam" "pVal_bam" "mCh_FPKM" "Sxl_FPKM" "bamF_FPKM" "bamM_FPKM" "Chr" "Start" "End" "Strand" "Attribute" "Transcript" > temp.txt
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "Gene_id" "Gene" "Motif" "FC_RNAi" "pVal_RNAi" "FC_bam" "pVal_bam" "mCh_FPKM" "Sxl_FPKM" "bamF_FPKM" "bamM_FPKM" "Chr" "Start" "End" "Strand" "Attribute" "Transcript" > temp1.txt
	FC=${i#/Users/cmdb/sxl_data/data_tables/${MOTIF}_gene_}
	FC=${FC%_annot_stats.txt}
	if [ $FC == "pos" ]
	then
		cat $i | sort -nrk4 | uniq >> temp1.txt
	else
		cat $i | sort -nk4 | uniq >> temp1.txt
	fi
	mv temp1.txt ${i%txt}sort_uniq.txt
	cat $i >> temp.txt; rm $i; mv temp.txt $i
done

rm ~/sxl_data/data_tables/${MOTIF}_exon_*not_support.txt
~/sxl_analysis/scripts/exon_annotation_compare.py ~/sxl_data/data_tables/${MOTIF}_exon_pos_annot_stats.txt annotate > ~/sxl_data/data_tables/${MOTIF}_exon_pos_annot_support.txt
~/sxl_analysis/scripts/exon_annotation_compare.py ~/sxl_data/data_tables/${MOTIF}_exon_neg_annot_stats.txt annotate > ~/sxl_data/data_tables/${MOTIF}_exon_neg_annot_support.txt
~/sxl_analysis/scripts/exon_annotation_compare.py ~/sxl_data/data_tables/${MOTIF}_exon_pos_annot_stats.txt unannotate > ~/sxl_data/data_tables/${MOTIF}_exon_pos_unannot_support.txt
~/sxl_analysis/scripts/exon_annotation_compare.py ~/sxl_data/data_tables/${MOTIF}_exon_neg_annot_stats.txt unannotate > ~/sxl_data/data_tables/${MOTIF}_exon_neg_unannot_support.txt
for i in ~/sxl_data/data_tables/${MOTIF}_exon_*not_support.txt
do
	FILE=${i#/Users/cmdb/sxl_data/data_tables/${MOTIF}_exon_}
	FILE=${FILE%_*nnot_support.txt}
	awk 'NR==1 {print}' $i > temp.txt
	if [ $FILE == "pos" ]
	then
		awk '{if ($20~"Yes") print}' $i | sort -nrk4 | uniq >> temp.txt; rm $i; mv temp.txt $i
	else
		awk '{if ($20~"Yes") print}' $i | sort -nk4 | uniq >> temp.txt; rm $i; mv temp.txt $i
	fi
done


rm ~/sxl_data/data_tables/${MOTIF}_novel_exon_*_unannot_stat*.txt
~/sxl_analysis/scripts/target_novel_exon.py stats ~/sxl_data/data_tables/novel_exon_pos_fc_unannot.txt > ~/sxl_data/data_tables/${MOTIF}_novel_exon_pos_unannot_stats.txt
~/sxl_analysis/scripts/target_novel_exon.py stats ~/sxl_data/data_tables/novel_exon_neg_fc_unannot.txt > ~/sxl_data/data_tables/${MOTIF}_novel_exon_neg_unannot_stats.txt

for i in ~/sxl_data/data_tables/${MOTIF}_novel_exon_*_unannot_stats.txt
do
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "Gene_id" "Gene" "Motif" "FC_RNAi" "pVal_RNAi" "FC_bam" "pVal_bam" "mCh_FPKM" "Sxl_FPKM" "bamF_FPKM" "bamM_FPKM" "Chr" "Start" "End" "Strand" "Exon_start" "Exon_end" "Prime" "e_id" > temp.txt
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "Gene_id" "Gene" "Motif" "FC_RNAi" "pVal_RNAi" "FC_bam" "pVal_bam" "mCh_FPKM" "Sxl_FPKM" "bamF_FPKM" "bamM_FPKM" "Chr" "Start" "End" "Strand" "exon_start" "exon_end" "Prime" "e_id" > temp1.txt
	FC=${i#/Users/cmdb/sxl_data/data_tables/${MOTIF}_novel_exon_}
	FC=${FC%_unannot_stats.txt}
	if [ $FC == "pos" ]
	then
		cat $i | sort -nrk4 | uniq >> temp1.txt
	else
		cat $i | sort -nk4 | uniq >> temp1.txt
	fi
	mv temp1.txt ${i%txt}sort_uniq.txt
	cat $i >> temp.txt; rm $i; mv temp.txt $i
done

rm ~/sxl_data/data_tables/${MOTIF}_novel_exons_*_unannot_supported.txt
~/sxl_analysis/scripts/exon_annotation_compare.py ~/sxl_data/data_tables/${MOTIF}_novel_exon_pos_unannot_stats.txt unannotate > ~/sxl_data/data_tables/${MOTIF}_novel_exons_pos_unannot_supported.txt
~/sxl_analysis/scripts/exon_annotation_compare.py ~/sxl_data/data_tables/${MOTIF}_novel_exon_neg_unannot_stats.txt unannotate > ~/sxl_data/data_tables/${MOTIF}_novel_exons_neg_unannot_supported.txt

for i in ~/sxl_data/data_tables/${MOTIF}_novel_exons_*_unannot_supported.txt
do
	FILE=${i#/Users/cmdb/sxl_data/data_tables/${MOTIF}_novel_exons_}
	FILE=${FILE%_unannot_supported.txt}
	awk 'NR==1 {print}' $i > temp.txt
	if [ $FILE == "pos" ]
	then
		awk '{if ($20~"Yes") print}' $i | sort -nrk4 | uniq >> temp.txt; rm $i; mv temp.txt $i
	else
		awk '{if ($20~"Yes") print}' $i | sort -nk4 | uniq >> temp.txt; rm $i; mv temp.txt $i
	fi
done
