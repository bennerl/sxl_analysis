#!/bin/bash
# 1 is exon file
# 2 is number of sequences to search for a motif
rm -r ~/meme_out
FILE=${1#/Users/cmdb/sxl_data/data_tables/exon_}
FILE=${FILE%_fc_annot.txt}

if [ $3 == "random" ]
then
	grep -vE "rDNA|mitochond|mapped_Scaff|2110000222" $1 > random.txt
	NUM=$(wc -l random.txt | sed "s|random.txt||g")
	python meme_random.py `echo $NUM` $2 random.txt > meme_random.txt
	/usr/local/opt/meme/bin/meme-chip -norc -oc ~/meme_out -db /Users/cmdb/annotations/motif_databases/JASPAR/JASPAR_CORE_2016_insects.meme meme_random.txt
	rm random.txt
	rm meme_random.txt
else
	if [ $FILE == "pos" ]
	then
		awk 'NR>1 {print}' $1 | sort -nrk12 | grep -vE "rDNA|mitochond|mapped_Scaff|2110000222" | head -${2} > temp.txt
	else
		awk 'NR>1 {print}' $1 | sort -nk12 | grep -vE "rDNA|mitochond|mapped_Scaff|2110000222" | head -${2} > temp.txt
	fi
	python meme.py > meme.txt
	/usr/local/opt/meme/bin/meme-chip -norc -oc ~/meme_out -db /Users/cmdb/annotations/motif_databases/JASPAR/JASPAR_CORE_2016_insects.meme meme.txt
	rm temp.txt
	rm meme.txt
fi


