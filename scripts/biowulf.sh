##########################################################
DIR=/data/bennerl/bam_sxl_RNAseq_fastq/
rm ${DIR}hisat_swarm.txt
for i in ${DIR}samples/*/
do
	FILE=${i#${DIR}samples/}
	FILE=${FILE%/}
	echo "hisat2 -x /data/bennerl/reference/r6.12_dmel_noERCC --fr --dta --known-splicesite-infile /data/bennerl/reference/r6.12_dmel_noERCC.splicesites.txt --novel-splicesite-outfile ${i}${FILE}_hisat_novel_splicesites.txt -1 ${i}*1.fastq -2 ${i}*2.fastq -p 8 -S ${i}${FILE}.sam" >> ${DIR}hisat_swarm.txt
done
swarm -g 16 -t 8 --time 01:00:00 --module hisat/2.0.4 -f ${DIR}hisat_swarm.txt
##########################################################
DIR=/data/bennerl/bam_sxl_RNAseq_fastq/
rm ${DIR}samtools_swarm.txt
for i in ${DIR}samples/*/
do
	FILE=${i#${DIR}samples/}
	FILE=${FILE%/}
	echo "samtools sort -@ 8 -o ${i}${FILE}.sorted.bam ${i}${FILE}.sam; samtools index ${i}${FILE}.sorted.bam; samtools sort -@ 8 -n -o ${i}${FILE}.name.sorted.bam ${i}${FILE}.sam" >> ${DIR}samtools_swarm.txt
done
swarm -g 16 -t 8 --time 01:00:00 --module samtools/1.3.1 -f ${DIR}samtools_swarm.txt
##########################################################
DIR=/data/bennerl/bam_sxl_RNAseq_fastq/
rm ${DIR}stringtie_gtf_swarm.txt
for i in ${DIR}samples/*/
do
	FILE=${i#${DIR}samples/}
	FILE=${FILE%/}
	echo "stringtie ${i}${FILE}.sorted.bam -o ${i}${FILE}.gtf -p 8 -G /data/bennerl/reference/r6.12_dmel_noERCC.gtf " >> ${DIR}stringtie_gtf_swarm.txt
done
swarm -g 16 -t 8 --time 00:30:00 --module stringtie/1.2.3 -f ${DIR}stringtie_gtf_swarm.txt
##########################################################
DIR=/data/bennerl/bam_sxl_RNAseq_fastq/
rm ${DIR}stringtie_merge_swarm.txt
rm ${DIR}cuffmerge_merge_swarm.txt
rm -r ${DIR}annotation/*
mkdir ${DIR}annotation/cuffmerge
mkdir ${DIR}annotation/stringtie
for i in ${DIR}samples/*/
do
	FILE=${i#${DIR}samples/}
	FILE=${FILE%/}
	cp ${i}${FILE}.gtf ${DIR}annotation/
	echo "${DIR}annotation/${FILE}.gtf" >> ${DIR}annotation/cuffmerge/cuffmerge_samples.txt
done
echo "stringtie --merge -G /data/bennerl/reference/r6.12_dmel_noERCC.gtf -o ${DIR}annotation/stringtie/stringtie_merge.gtf ${DIR}annotation/*gtf" >> ${DIR}stringtie_merge_swarm.txt
echo "cuffmerge -g /data/bennerl/reference/r6.12_dmel_noERCC.gtf -p 8 -s /data/bennerl/reference/r6.12_dmel_noERCC.fasta -o ${DIR}annotation/cuffmerge ${DIR}annotation/cuffmerge/cuffmerge_samples.txt" >> ${DIR}cuffmerge_merge_swarm.txt
swarm -g 16 -t 8 --time 00:30:00 --module cufflinks/2.2.1_patched -f ${DIR}cuffmerge_merge_swarm.txt
swarm -g 16 -t 8 --time 00:30:00 --module stringtie/1.2.3 -f ${DIR}stringtie_merge_swarm.txt
##########################################################
DIR=/data/bennerl/bam_sxl_RNAseq_fastq/
rm ${DIR}stringtie_ballgown_swarm.txt
for i in ${DIR}samples/*/
do
	rm -r ${i}stringtie_annotated
	rm -r ${i}stringtie_merge
	rm -r ${i}stringtie_cuffmerge
	FILE=${i#${DIR}samples/}
	FILE=${FILE%/}
	mkdir ${i}stringtie_annotated
	mkdir ${i}stringtie_merge
	mkdir ${i}stringtie_cuffmerge
	echo "stringtie ${i}${FILE}.sorted.bam -o ${i}/stringtie_annotated/${FILE}.gtf -A ${i}/stringtie_annotated/${FILE}_gene_abundances.tab -B -e -p 8 -G /data/bennerl/reference/r6.12_dmel_noERCC.gtf " >> ${DIR}stringtie_ballgown_swarm.txt
	echo "stringtie ${i}${FILE}.sorted.bam -o ${i}/stringtie_merge/${FILE}.gtf -A ${i}/stringtie_merge/${FILE}_gene_abundances.tab -B -e -p 8 -G ${DIR}annotation/stringtie/stringtie_merge.gtf" >> ${DIR}stringtie_ballgown_swarm.txt
	echo "stringtie ${i}${FILE}.sorted.bam -o ${i}/stringtie_cuffmerge/${FILE}.gtf -A ${i}/stringtie_cuffmerge/${FILE}_gene_abundances.tab -B -e -p 8 -G ${DIR}annotation/cuffmerge/merged.gtf" >> ${DIR}stringtie_ballgown_swarm.txt
done
swarm -g 16 -t 8 --time 00:30:00 --module stringtie/1.2.3 -f ${DIR}stringtie_ballgown_swarm.txt
##########################################################
rm -r ~/ballgown
mkdir ~/ballgown
DIR=/data/bennerl/bam_sxl_RNAseq_fastq/
mkdir ~/ballgown/stringtie_annotated
mkdir ~/ballgown/stringtie_merge
mkdir ~/ballgown/stringtie_cuffmerge
for i in ${DIR}samples/*/
do
	FILE=${i#${DIR}samples/}
	FILE=${FILE%/}
	mkdir ~/ballgown/stringtie_annotated/${FILE}
	mkdir ~/ballgown/stringtie_merge/${FILE}
	mkdir ~/ballgown/stringtie_cuffmerge/${FILE}
	cp ${i}/stringtie_annotated/*ctab ~/ballgown/stringtie_annotated/${FILE}
	cp ${i}/stringtie_merge/*ctab ~/ballgown/stringtie_merge/${FILE}
	cp ${i}/stringtie_cuffmerge/*ctab ~/ballgown/stringtie_cuffmerge/${FILE}
done
##########################################################











