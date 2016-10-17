#!/bin/bash
#for running on biowulf, paired end only

#a lot of the for loops are redundant file checks with standard software exec commands, prob best to clean up by defining one function to handle all these tasks, not sure exactly howd that be done though because each loop has slightly different variables and setting, might be good to just keep it modular like this

#need to alter annotation indexes wtih ERCCs for hisat2 according to stringtie pipeline, i dont think it is optimized for determining novel splicesites, does not have a knownsplicesites file to read into.

#also need to doublecheck hisat2 command according to stringtie pipeline, should outline good pipeline according to paper, ie determine novel splicesites and then remap with these novel splicesites (i think)

#i believe i have the option of creating a robust novel splicesite file according to treatment.  If i --merge novel splicesites found in individual smaples with stringtie according to genotype, then it should give me an output of these that are

#WORKFLOW TO ADJUST
# hisat align with splicesites -> run stringtie with -G -> stringtie --merge with

#look further into running cuffmerge and identifying novel transcripts that are found or not found betwwn treatments, specifically sxl- vs sxl+ and bam- males and bam- females

#generate QC data for samples

#use gffcompare to report entirely novel

#you can re-run hisat2 and provide it with the novel splicsites you generated in the first run, does this increase mapping for downstream application like stringtie? or does this not matter for stringtie?
###############################################################
#first argument should be the path of the fastq folders containing separate fastq files
if [ -z $1 ]
then
	echo ""
	echo "You didnt supply the path to the samples"
	echo ""
	exit
fi
###############################################################
#set absolute path of fastq folders
DIR=$1
EXEC_HOME=`pwd`
###############################################################
#check to see if fasta file is compressed, if not, it will uncompress, needs to be in .gz format obviously
if [ -f ${DIR}gunzip_swarm.txt ]
then
	rm ${DIR}gunzip_swarm.txt
fi

for i in ${DIR}*/*.fast*
do
	if [ ! -f ${i%.gz} ]
	then
		echo "gunzip $i" >> ${DIR}gunzip_swarm.txt
	else
		echo "${i%.gz} exists"
	fi
done

if [ -f ${DIR}gunzip_swarm.txt ]
then
	JOB_GUNZIP=$(swarm -g 32 -t 8 --time 00:30:00 -f ${DIR}gunzip_swarm.txt)
	echo $JOB_GUNZIP
else
	echo "Not running ${DIR}gunzip_swarm.txt"
fi
###############################################################
#module load hisat/2.0.4
if [ -f ${DIR}hisat_swarm.txt ]
then
	rm ${DIR}hisat_swarm.txt
fi

for i in ${DIR}*/
do
	FILE=${i#${DIR}}
	FILE=${FILE%/}
	if [ ! -f ${i}${FILE}.sam ]
	then
		MATE_1S=`echo ${i}*_1.fastq`
		MATE_2S=`echo ${i}*_2.fastq`
		echo "hisat2 -x /data/bennerl/reference/r6.12_dmel_noERCC --fr --dta --known-splicesite-infile /data/bennerl/reference/r6.12_dmel_noERCC.splicesites.txt --novel-splicesite-outfile ${i}${FILE}_hisat_novel_splicesites.txt -1 $MATE_1S -2 $MATE_2S -p 8 -S ${i}${FILE}.sam" >> ${DIR}hisat_swarm.txt
	else
		echo "File ${i}${FILE}.sam exists"
	fi
done

if [ -f ${DIR}hisat_swarm.txt ]
then
	if [ ! -z $JOB_GUNZIP ]
	then
		JOB_HISAT=$(swarm -g 32 -t 8 --time 01:00:00 --module hisat/2.0.4 -f ${DIR}hisat_swarm.txt --dependency afterany:$JOB_GUNZIP)
		echo $JOB_HISAT
	else
		JOB_HISAT=$(swarm -g 32 -t 8 --time 01:00:00 --module hisat/2.0.4 -f ${DIR}hisat_swarm.txt)
		echo $JOB_HISAT
	fi
else
	echo "Not running ${DIR}hisat_swarm.txt"
fi
###############################################################
#module load samtools/1.3.1
if [ -f ${DIR}samtools_swarm.txt ]
then
	rm ${DIR}samtools_swarm.txt
fi

for i in ${DIR}*/
do
	FILE=${i#${DIR}}
	FILE=${FILE%/}
	if [ ! -f ${i}${FILE}.sorted.bam ]
	then
		echo "samtools sort -o ${i}${FILE}.sorted.bam ${i}${FILE}.sam; samtools index ${i}${FILE}.sorted.bam; samtools sort -n -o ${i}${FILE}.name.sorted.bam ${i}${FILE}.sam" >> ${DIR}samtools_swarm.txt
	else
		echo "File ${i}${FILE}.sorted.bam exists"
	fi
done

if [ -f ${DIR}samtools_swarm.txt ]
then
	if [ ! -z $JOB_HISAT ]
	then
		JOB_SAMTOOL=$(swarm -g 32 -t 8 --time 01:00:00 --module samtools/1.3.1 -f ${DIR}samtools_swarm.txt --dependency afterany:$JOB_HISAT)
		echo $JOB_SAMTOOL
	else
		JOB_SAMTOOL=$(swarm -g 32 -t 8 --time 01:00:00 --module samtools/1.3.1 -f ${DIR}samtools_swarm.txt)
		echo $JOB_SAMTOOL
	fi
else
	echo "Not running ${DIR}samtools_swarm.txt"
fi
###############################################################
#module load stringtie/1.2.3
#this is the first stringtie command generating novel_spicesite GTFs to use for merge
if [ -f ${DIR}stringtie_GTF_swarm.txt ]
then
	rm ${DIR}stringtie_GTF_swarm.txt
fi

for i in ${DIR}*/
do
	FILE=${i#${DIR}}
	FILE=${FILE%/}
	if [ ! -f ${i}${FILE}.gtf ]
	then
		echo "stringtie ${i}${FILE}.sorted.bam -o ${i}${FILE}.gtf -p 8 -G /data/bennerl/reference/r6.12_dmel_noERCC.gtf " >> ${DIR}stringtie_GTF_swarm.txt
	else
		echo "File ${i}${FILE}.gtf exists"
	fi
done

if [ -f ${DIR}stringtie_GTF_swarm.txt ]
then
	if [ ! -z $JOB_SAMTOOL ]
	then
		JOB_STRINGTIE_GTF=$(swarm -g 32 -t 8 --time 00:45:00 --module stringtie/1.2.3 -f ${DIR}stringtie_GTF_swarm.txt --dependency afterany:$JOB_SAMTOOL)
		echo $JOB_STRINGTIE_GTF
	else
		JOB_STRINGTIE_GTF=$(swarm -g 32 -t 8 --time 00:45:00 --module stringtie/1.2.3 -f ${DIR}stringtie_GTF_swarm.txt)
		echo $JOB_STRINGTIE_GTF
	fi
else
	echo "Not running ${DIR}stringtie_GTF_swarm.txt"
fi
###############################################################
#stringtie --merge

TREATMENTS=()
for i in ${DIR}*/
do
	FILE=${i#${DIR}}
	FILE=${FILE%/}
	TREATMENTS+=(${FILE%-rep*})
done
TREATMENTS=$(echo ${TREATMENTS[@]} | tr ' ' '\n' | sort | uniq)

if [ ! -d ${DIR}../annotation/stringtie_all_sample_annotation/ ]
then
	mkdir ${DIR}../annotation/stringtie_all_sample_annotation/
fi

for TREAT in $TREATMENTS
do
	if [ ! -d ${DIR}../annotation/${TREAT}_merge/ ]
	then
		mkdir ${DIR}../annotation/${TREAT}_merge/
	fi
	for i in ${DIR}*/
	do
		FILE=${i#${DIR}}
		FILE=${FILE%/}
		if [ ${FILE%-rep*} == $TREAT ]
		then
			if [ -f ${i}${FILE}.gtf ]
			then
				cp ${i}${FILE}.gtf ${DIR}../annotation/${TREAT}_merge/
				cp ${i}${FILE}.gtf ${DIR}../annotation/stringtie_all_sample_annotation/
			fi
		fi
	done
done

if [ -f ${DIR}stringtie_merge_swarm.txt ]
then
	rm ${DIR}stringtie_merge_swarm.txt
fi

for TREAT in ${DIR}../annotation/*_merge/
do
	FILE=${TREAT#${DIR}../annotation/}
	FILE=${FILE%_merge/}
	if [ ! -f ${TREAT}${FILE}_merge.gtf ]
	then
		echo "stringtie --merge -G /data/bennerl/reference/r6.12_dmel_noERCC.gtf -o ${TREAT}${FILE}_merge.gtf ${TREAT}${FILE}*gtf" >> ${DIR}stringtie_merge_swarm.txt
	fi
done

if [ ! -f ${DIR}../annotation/stringtie_all_sample_annotation/new_annotation_merge.gtf ]
then
	echo "stringtie --merge -G /data/bennerl/reference/r6.12_dmel_noERCC.gtf -o ${DIR}../annotation/stringtie_all_sample_annotation/new_annotation_merge.gtf ${DIR}../annotation/stringtie_all_sample_annotation/*gtf" >> ${DIR}stringtie_merge_swarm.txt
fi

if [ -f ${DIR}stringtie_merge_swarm.txt ]
then
	if [ ! -z $JOB_STRINGTIE_GTF ]
	then
		JOB_STRINGTIE_MERGE=$(swarm -g 32 -t 8 --time 01:00:00 --module stringtie/1.2.3 -f ${DIR}stringtie_merge_swarm.txt --dependency afterany:$JOB_STRINGTIE_GTF)
		echo $JOB_STRINGTIE_MERGE
	else
		JOB_STRINGTIE_MERGE=$(swarm -g 32 -t 8 --time 01:00:00 --module stringtie/1.2.3 -f ${DIR}stringtie_merge_swarm.txt)
		echo $JOB_STRINGTIE_MERGE
	fi
else
	echo "Not running ${DIR}stringtie_merge_swarm.txt"
fi

##############################################################
#need to add cuffcompare command to merge gtf files made in different treatments so that all transcripts are found in every alignement file

if [ ! -d ${DIR}../annotation/gffcompare/ ]
then
	mkdir ${DIR}../annotation/gffcompare
	mkdir ${DIR}../annotation/gffcompare/merge_all_treat
fi
f
if [ -f ${DIR}gffcompare_all_treat_swarm.txt ]
then
	rm ${DIR}gffcompare_all_treat_swarm.txt
fi

if [ ! -f ${DIR}../annotation/gffcompare/merge_all_treat/merged_all_treat.combined.gtf ]
then
	for i in ${DIR}../annotation/*_merge/
	do
		#this errors because I need to write in a filecheck, but for now its fine
		cp ${i}*merge.gtf ${i}../gffcompare/merge_all_treat/
	done
	if [ ! -z $JOB_STRINGTIE_MERGE ]
	then
		echo "gffcompare -r /data/Oliverlab/data/FlyBase/FB2015_04/dmel_ERCC.gtf -o ${DIR}../annotation/gffcompare/merge_all_treat/merged_all_treat ${DIR}../annotation/gffcompare/merge_all_treat/*.gtf" > ${DIR}gffcompare_all_treat_swarm.txt
		JOB_GFF_COMPARE_ALL_TREAT=$(swarm -g 32 -t 8 --time 01:00:00 -f ${DIR}gffcompare_all_treat_swarm.txt --dependency afterany:$JOB_STRINGTIE_MERGE)
	else
		echo "gffcompare -r /data/Oliverlab/data/FlyBase/FB2015_04/dmel_ERCC.gtf -o ${DIR}../annotation/gffcompare/merge_all_treat/merged_all_treat ${DIR}../annotation/gffcompare/merge_all_treat/*.gtf" > ${DIR}gffcompare_all_treat_swarm.txt
		JOB_GFF_COMPARE_ALL_TREAT=$(swarm -g 32 -t 8 --time 01:00:00 -f ${DIR}gffcompare_all_treat_swarm.txt)
	fi
fi
###############################################################

#strintie and ballgown

#module load stringtie/1.2.3
#this is the second stringtie command generating novel_spicesite GTFs to use for merge
if [ -f ${DIR}stringtie_ballgown_swarm.txt ]
then
	rm ${DIR}stringtie_ballgown_swarm.txt
fi

for i in ${DIR}*/
do
	FILE=${i#${DIR}}
	FILE=${FILE%/}
	if [ ! -d ${i}ballgown/ ]
	then
		mkdir ${i}ballgown/
	fi
	if [ ! -f ${i}ballgown/${FILE}.gtf ]
	then
		echo "stringtie ${i}${FILE}.sorted.bam -o ${i}/ballgown/${FILE}.gtf -A ${i}/ballgown/${FILE}_gene_abundances.tab -B -e -p 8 -G ${DIR}../annotation/stringtie_all_sample_annotation/new_annotation_merge.gtf " >> ${DIR}stringtie_ballgown_swarm.txt
	else
		echo "File ${i}ballgown/${FILE}.gtf exists"
	fi
done

if [ -f ${DIR}stringtie_ballgown_swarm.txt ]
then
	if [ ! -z $JOB_STRINGTIE_MERGE ]
	then
		JOB_STRINGTIE_BALLGOWN=$(swarm -g 32 -t 8 --time 00:30:00 --module stringtie/1.2.3 -f ${DIR}stringtie_ballgown_swarm.txt --dependency afterany:$JOB_STRINGTIE_MERGE)
		echo $JOB_STRINGTIE_BALLGOWN
	else
		JOB_STRINGTIE_BALLGOWN=$(swarm -g 32 -t 8 --time 00:30:00 --module stringtie/1.2.3 -f ${DIR}stringtie_ballgown_swarm.txt)
		echo $JOB_STRINGTIE_BALLGOWN
	fi
else
	echo "Not running ${DIR}stringtie_ballgown_swarm.txt"
fi

###############################################################
#htseq-count
DIR=$(pwd)
for i in /data/bennerl/bam_sxl_RNAseq_fastq/samples/*/
do
	cd $i
	FILE=${i#/data/bennerl/bam_sxl_RNAseq_fastq/samples/}
	FILE=${FILE%/}
	echo "htseq-count -f bam -r pos -s yes ${i}*sorted.bam /data/bennerl/reference/r6.12_dmel_noERCC.gtf > ${i}${FILE}_htseq-counts_posY.txt" >> /data/bennerl/bam_sxl_RNAseq_fastq/samples/htseq-count_swarm.txt
	cd $DIR
done

swarm -g 32 -t 8 --time 02:00:00 --module htseq/0.6.1p1 -f /data/bennerl/bam_sxl_RNAseq_fastq/samples/htseq-count_swarm.txt


##############################################################
#fastqc
module load fastqc/0.11.5

for i in /data/bennerl/bam_sxl_RNAseq_fastq/samples/*/
do
	cd $i
	fastqc -o /data/bennerl/bam_sxl_RNAseq_fastq/fastqc_reports/ *fastq
	cd
done
##############################################################

echo "stringtie ${i}${FILE}.sorted.bam -A ${i}${FILE}.annotated_abundances.tab -b ${i} -e -p 8 -G /data/bennerl/reference/r6.12_dmel_noERCC.gtf " >> /data/bennerl/bam_sxl_RNAseq_fastq/samples/stringtie_GTF_swarm.txt

##############################################################

for i in /data/bennerl/bam_sxl_RNAseq_fastq/samples/*/
do
	FILE=${i#/data/bennerl/bam_sxl_RNAseq_fastq/samples/}; FILE=${FILE%/}
	cp ${i}${FILE}.gtf /data/bennerl/bam_sxl_RNAseq_fastq/annotation
done
echo "stringtie --merge -G /data/bennerl/reference/r6.12_dmel_noERCC.gtf -o /data/bennerl/bam_sxl_RNAseq_fastq/annotation/bam-RNAi_new_annotation_merge.gtf /data/bennerl/bam_sxl_RNAseq_fastq/annotation/*gtf" > /data/bennerl/bam_sxl_RNAseq_fastq/samples/stringtie_merge_swarm.txt
swarm -g 32 -t 8 --time 00:30:00 --module stringtie/1.2.3 -f /data/bennerl/bam_sxl_RNAseq_fastq/samples/stringtie_merge_swarm.txt
##############################################################

echo "cuffmerge -g /data/bennerl/reference/r6.12_dmel_noERCC.gtf -p 8 -s /data/bennerl/reference/r6.12_dmel_noERCC.fasta -o /data/bennerl/bam_sxl_RNAseq_fastq/annotation/ /data/bennerl/bam_sxl_RNAseq_fastq/samples/cuffmerge_samples.txt" > /data/bennerl/bam_sxl_RNAseq_fastq/samples/cuffmerge_swarm.txt
swarm -g 32 -t 8 --time 00:30:00 --module cufflinks/2.2.1_patched -f /data/bennerl/bam_sxl_RNAseq_fastq/samples/cuffmerge_swarm.txt

##############################################################

for i in /data/bennerl/bam_sxl_RNAseq_fastq/samples/*/
do
	FILE=${i#/data/bennerl/bam_sxl_RNAseq_fastq/samples/}; FILE=${FILE%/}
	mkdir ${i}merged_ballgown_stringtie
	mkdir ${i}merged_ballgown_cuffmerge
	echo "stringtie ${i}${FILE}.sorted.bam -o ${i}/merged_ballgown_stringtie/${FILE}.merged.gtf -A ${i}/merged_ballgown_stringtie/${FILE}_merged_gene_abundances.tab -B -e -p 8 -G /data/bennerl/bam_sxl_RNAseq_fastq/annotation/bam-RNAi_new_annotation_merge.gtf " >> /data/bennerl/bam_sxl_RNAseq_fastq/samples/stringtie_ballgown_swarm.txt
	echo "stringtie ${i}${FILE}.sorted.bam -o ${i}/merged_ballgown_cuffmerge/${FILE}.merged.gtf -A ${i}/merged_ballgown_cuffmerge/${FILE}_merged_gene_abundances.tab -B -e -p 8 -G /data/bennerl/bam_sxl_RNAseq_fastq/annotation/merge.gtf " >> /data/bennerl/bam_sxl_RNAseq_fastq/samples/stringtie_ballgown_swarm.txt
done
swarm -g 32 -t 8 --time 00:30:00 --module stringtie/1.2.3 -f /data/bennerl/bam_sxl_RNAseq_fastq/samples/stringtie_ballgown_swarm.txt
##############################################################

for i in /data/bennerl/bam_sxl_RNAseq_fastq/samples/*/; do FILE=${i#/data/bennerl/bam_sxl_RNAseq_fastq/samples/}; FILE=${FILE%/}; mkdir ~/ballgown/merged_stringtie/${FILE}; cp ${i}/merged_ballgown_stringtie/*ctab ~/ballgown/merged_stringtie/${FILE}; done
