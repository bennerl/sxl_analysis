#change motif in target_finder.py and 
MOTIF="Sakashita"
RScript /Users/cmdb/bam_sxl_analysis/scripts/exon_gene_select.R #
RScript /Users/cmdb/bam_sxl_analysis/scripts/exon_unannotated_gene_select.R #

rm /Users/cmdb/bam_sxl_analysis/data_tables/gene_*_annot.chr.txt
for i in /Users/cmdb/bam_sxl_analysis/data_tables/gene_*_annot.txt
do
	/Users/cmdb/bam_sxl_analysis/scripts/add_chr_to_BG_gene_list.py $i > temp.txt
	mv temp.txt ${i%txt}chr.txt
done

rm ~/bam_sxl_analysis/bed/${MOTIF}_exon_*_annot.bed
~/bam_sxl_analysis/scripts/target_finder.py exon bed ~/bam_sxl_analysis/data_tables/exon_pos_fc_annot.txt ~/bam_sxl_analysis/data_tables/e2g_annotated.txt > ~/bam_sxl_analysis/bed/${MOTIF}_exon_pos_annot.bed
~/bam_sxl_analysis/scripts/target_finder.py exon bed ~/bam_sxl_analysis/data_tables/exon_neg_fc_annot.txt ~/bam_sxl_analysis/data_tables/e2g_annotated.txt > ~/bam_sxl_analysis/bed/${MOTIF}_exon_neg_annot.bed

for i in ~/bam_sxl_analysis/bed/${MOTIF}_exon_*_annot.bed
do
	NAME=${i#/Users/cmdb/bam_sxl_analysis/bed/}
	NAME=${NAME%.bed}
	NAME_SPACE=$(echo $NAME | tr '_' ' ')
	NAME=$(echo $NAME | sed 's|_||g')
	echo "track name=\"$NAME_SPACE\" description=\"$NAME\" visibility=2 itemRgb=\"On\"" > temp.txt
	cat $i | sort | uniq >> temp.txt; rm $i; mv temp.txt $i
done

rm ~/bam_sxl_analysis/bed/${MOTIF}_gene*annot.bed
~/bam_sxl_analysis/scripts/target_finder.py gene bed ~/bam_sxl_analysis/data_tables/gene_pos_fc_annot.chr.txt ~/annotations/r6.12_dmel_noERCC.gtf > ~/bam_sxl_analysis/bed/${MOTIF}_gene_pos_annot.bed
~/bam_sxl_analysis/scripts/target_finder.py gene bed ~/bam_sxl_analysis/data_tables/gene_neg_fc_annot.chr.txt ~/annotations/r6.12_dmel_noERCC.gtf > ~/bam_sxl_analysis/bed/${MOTIF}_gene_neg_annot.bed

for i in ~/bam_sxl_analysis/bed/${MOTIF}_gene*annot.bed
do
	NAME=${i#/Users/cmdb/bam_sxl_analysis/bed/}
	NAME=${NAME%.bed}
	NAME_SPACE=$(echo $NAME | tr '_' ' ')
	NAME=$(echo $NAME | sed 's|_||g')
	echo "track name=\"$NAME_SPACE\" description=\"$NAME\" visibility=2 itemRgb=\"On\"" > temp.txt
	cat $i | sort | uniq >> temp.txt; rm $i; mv temp.txt $i
done

rm ~/bam_sxl_analysis/data_tables/novel_exons.txt
~/bam_sxl_analysis/scripts/novel_exon_find.py gtf > ~/bam_sxl_analysis/data_tables/novel_exons.txt
rm ~/bam_sxl_analysis/data_tables/novel_exon_*_fc_unannot.txt
~/bam_sxl_analysis/scripts/novel_exon_find.py merge ~/bam_sxl_analysis/data_tables/exon_pos_fc_unannot.txt > ~/bam_sxl_analysis/data_tables/novel_exon_pos_fc_unannot.txt
~/bam_sxl_analysis/scripts/novel_exon_find.py merge ~/bam_sxl_analysis/data_tables/exon_neg_fc_unannot.txt > ~/bam_sxl_analysis/data_tables/novel_exon_neg_fc_unannot.txt
rm ~/bam_sxl_analysis/bed/${MOTIF}_novel_exon_*_unannot.bed
~/bam_sxl_analysis/scripts/target_novel_exon.py bed ~/bam_sxl_analysis/data_tables/novel_exon_pos_fc_unannot.txt > ~/bam_sxl_analysis/bed/${MOTIF}_novel_exon_pos_unannot.bed
~/bam_sxl_analysis/scripts/target_novel_exon.py bed ~/bam_sxl_analysis/data_tables/novel_exon_neg_fc_unannot.txt > ~/bam_sxl_analysis/bed/${MOTIF}_novel_exon_neg_unannot.bed

for i in ~/bam_sxl_analysis/bed/${MOTIF}_novel*.bed
do
	NAME=${i#/Users/cmdb/bam_sxl_analysis/bed/}
	NAME=${NAME%.bed}
	NAME_SPACE=$(echo $NAME | tr '_' ' ')
	NAME=$(echo $NAME | sed 's|_||g')
	echo "track name=\"$NAME_SPACE\" description=\"$NAME\" visibility=2 itemRgb=\"On\"" > temp.txt
	cat $i | sort | uniq >> temp.txt; rm $i; mv temp.txt $i
done

rm ~/bam_sxl_analysis/data_tables/${MOTIF}_exon_*_annot_stat*.txt
~/bam_sxl_analysis/scripts/target_finder.py exon stats ~/bam_sxl_analysis/data_tables/exon_pos_fc_annot.txt ~/bam_sxl_analysis/data_tables/e2g_annotated.txt > ~/bam_sxl_analysis/data_tables/${MOTIF}_exon_pos_annot_stats.txt
~/bam_sxl_analysis/scripts/target_finder.py exon stats ~/bam_sxl_analysis/data_tables/exon_neg_fc_annot.txt ~/bam_sxl_analysis/data_tables/e2g_annotated.txt > ~/bam_sxl_analysis/data_tables/${MOTIF}_exon_neg_annot_stats.txt

for i in ~/bam_sxl_analysis/data_tables/${MOTIF}_exon_*stat*.txt
do
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "Gene_id" "Gene" "Motif" "FC_RNAi" "pVal_RNAi" "FC_bam" "pVal_bam" "mCh_FPKM" "Sxl_FPKM" "bamF_FPKM" "bamM_FPKM" "Chr" "Start" "End" "Strand" "Exon_start" "Exon_end" "Prime" "e_id" > temp.txt
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "Gene_id" "Gene" "Motif" "FC_RNAi" "pVal_RNAi" "FC_bam" "pVal_bam" "mCh_FPKM" "Sxl_FPKM" "bamF_FPKM" "bamM_FPKM" "Chr" "Start" "End" "Strand" "exon_start" "exon_end" "Prime" "e_id" > temp1.txt
	FC=${i#/Users/cmdb/bam_sxl_analysis/data_tables/${MOTIF}_exon_}
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

rm ~/bam_sxl_analysis/data_tables/${MOTIF}_gene_*_annot_stat*.txt
~/bam_sxl_analysis/scripts/target_finder.py gene stats ~/bam_sxl_analysis/data_tables/gene_pos_fc_annot.chr.txt ~/annotations/r6.12_dmel_noERCC.gtf > ~/bam_sxl_analysis/data_tables/${MOTIF}_gene_pos_annot_stats.txt
~/bam_sxl_analysis/scripts/target_finder.py gene stats ~/bam_sxl_analysis/data_tables/gene_neg_fc_annot.chr.txt ~/annotations/r6.12_dmel_noERCC.gtf > ~/bam_sxl_analysis/data_tables/${MOTIF}_gene_neg_annot_stats.txt

for i in ~/bam_sxl_analysis/data_tables/${MOTIF}_gene*stat*.txt
do
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "Gene_id" "Gene" "Motif" "FC_RNAi" "pVal_RNAi" "FC_bam" "pVal_bam" "mCh_FPKM" "Sxl_FPKM" "bamF_FPKM" "bamM_FPKM" "Chr" "Start" "End" "Strand" "Attribute" "Transcript" > temp.txt
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "Gene_id" "Gene" "Motif" "FC_RNAi" "pVal_RNAi" "FC_bam" "pVal_bam" "mCh_FPKM" "Sxl_FPKM" "bamF_FPKM" "bamM_FPKM" "Chr" "Start" "End" "Strand" "Attribute" "Transcript" > temp1.txt
	FC=${i#/Users/cmdb/bam_sxl_analysis/data_tables/${MOTIF}_gene_}
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

rm ~/bam_sxl_analysis/data_tables/${MOTIF}_exon_*not_support.txt
~/bam_sxl_analysis/scripts/exon_annotation_compare.py ~/bam_sxl_analysis/data_tables/${MOTIF}_exon_pos_annot_stats.txt annotate > ~/bam_sxl_analysis/data_tables/${MOTIF}_exon_pos_annot_support.txt
~/bam_sxl_analysis/scripts/exon_annotation_compare.py ~/bam_sxl_analysis/data_tables/${MOTIF}_exon_neg_annot_stats.txt annotate > ~/bam_sxl_analysis/data_tables/${MOTIF}_exon_neg_annot_support.txt
~/bam_sxl_analysis/scripts/exon_annotation_compare.py ~/bam_sxl_analysis/data_tables/${MOTIF}_exon_pos_annot_stats.txt unannotate > ~/bam_sxl_analysis/data_tables/${MOTIF}_exon_pos_unannot_support.txt
~/bam_sxl_analysis/scripts/exon_annotation_compare.py ~/bam_sxl_analysis/data_tables/${MOTIF}_exon_neg_annot_stats.txt unannotate > ~/bam_sxl_analysis/data_tables/${MOTIF}_exon_neg_unannot_support.txt
for i in ~/bam_sxl_analysis/data_tables/${MOTIF}_exon_*not_support.txt
do
	FILE=${i#/Users/cmdb/bam_sxl_analysis/data_tables/${MOTIF}_exon_}
	FILE=${FILE%_*nnot_support.txt}
	awk 'NR==1 {print}' $i > temp.txt
	if [ $FILE == "pos" ]
	then
		awk '{if ($20~"Yes") print}' $i | sort -nrk4 | uniq >> temp.txt; rm $i; mv temp.txt $i
	else
		awk '{if ($20~"Yes") print}' $i | sort -nk4 | uniq >> temp.txt; rm $i; mv temp.txt $i
	fi
done


rm ~/bam_sxl_analysis/data_tables/${MOTIF}_novel_exon_*_unannot_stat*.txt
~/bam_sxl_analysis/scripts/target_novel_exon.py stats ~/bam_sxl_analysis/data_tables/novel_exon_pos_fc_unannot.txt > ~/bam_sxl_analysis/data_tables/${MOTIF}_novel_exon_pos_unannot_stats.txt
~/bam_sxl_analysis/scripts/target_novel_exon.py stats ~/bam_sxl_analysis/data_tables/novel_exon_neg_fc_unannot.txt > ~/bam_sxl_analysis/data_tables/${MOTIF}_novel_exon_neg_unannot_stats.txt

for i in ~/bam_sxl_analysis/data_tables/${MOTIF}_novel_exon_*_unannot_stats.txt
do
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "Gene_id" "Gene" "Motif" "FC_RNAi" "pVal_RNAi" "FC_bam" "pVal_bam" "mCh_FPKM" "Sxl_FPKM" "bamF_FPKM" "bamM_FPKM" "Chr" "Start" "End" "Strand" "Exon_start" "Exon_end" "Prime" "e_id" > temp.txt
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "Gene_id" "Gene" "Motif" "FC_RNAi" "pVal_RNAi" "FC_bam" "pVal_bam" "mCh_FPKM" "Sxl_FPKM" "bamF_FPKM" "bamM_FPKM" "Chr" "Start" "End" "Strand" "exon_start" "exon_end" "Prime" "e_id" > temp1.txt
	FC=${i#/Users/cmdb/bam_sxl_analysis/data_tables/${MOTIF}_novel_exon_}
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

rm ~/bam_sxl_analysis/data_tables/${MOTIF}_novel_exons_*_unannot_supported.txt
~/bam_sxl_analysis/scripts/exon_annotation_compare.py ~/bam_sxl_analysis/data_tables/${MOTIF}_novel_exon_pos_unannot_stats.txt unannotate > ~/bam_sxl_analysis/data_tables/${MOTIF}_novel_exons_pos_unannot_supported.txt
~/bam_sxl_analysis/scripts/exon_annotation_compare.py ~/bam_sxl_analysis/data_tables/${MOTIF}_novel_exon_neg_unannot_stats.txt unannotate > ~/bam_sxl_analysis/data_tables/${MOTIF}_novel_exons_neg_unannot_supported.txt

for i in ~/bam_sxl_analysis/data_tables/${MOTIF}_novel_exons_*_unannot_supported.txt
do
	FILE=${i#/Users/cmdb/bam_sxl_analysis/data_tables/${MOTIF}_novel_exons_}
	FILE=${FILE%_unannot_supported.txt}
	awk 'NR==1 {print}' $i > temp.txt
	if [ $FILE == "pos" ]
	then
		awk '{if ($20~"Yes") print}' $i | sort -nrk4 | uniq >> temp.txt; rm $i; mv temp.txt $i
	else
		awk '{if ($20~"Yes") print}' $i | sort -nk4 | uniq >> temp.txt; rm $i; mv temp.txt $i
	fi
done
