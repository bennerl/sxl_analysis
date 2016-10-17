DIR=$(pwd)
cd /Users/cmdb/bam_ballgown/
for i in Sxl*
do
	mv $i ${i%r*}RNAi${i#*l}
done

for i in mCh*
do
	mv $i ${i%r*}RNAi${i#*h}
done
cd $DIR