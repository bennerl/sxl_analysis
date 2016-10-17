#Workflow
The workflow that I've completed so far is outlined below, actual commands are in the scripts/master.sh file.

For all sample, abbreviations are as follows
BamF/F -> Bam- females
BamM/M -> Bam- males
mCh/C -> Bam- females nanos-gal4 > mCh-RNAi
Sxl/S -> Bam- females nanos-gal4 > Sxl-RNAi

Map reads with Hisat2 ->
Generate guided splicing/gene annotation alignement files (GTF) for each sample with stringtie ->
Merge each samples GTF file with stringtie merge to generate an annotation file with characterized and novel gene/splicesites ->
Generate read alignements and transcript/gene abundance files with stringtie to use with ballgown ->
Determine differential gene expression and transcript abundance with ballgown.

#Plots
All plots/figures follow the same sample abbreviations as above

#MinIon
Total RNA from above samples were bioanalyzed (need to ask Caitlin for exact machine/protocl used) according to standard protocols (plots/total_RNA_bioanalyzer.pdf).  rRNA bands look great in RNAi samples, a little degraded in bam samples, but not that bad (i've sequenced much worse).  Measured RNA conc. with Qubit HS kit. Sample conc. are as follows:

Sample     -> Conc.      | uL pooled
Sxl-RNAi-1 -> 43.2 ng/uL | 11.6
Sxl-RNAi-2 -> 36.6 ng/uL | 13.7
Sxl-RNAi-3 -> 45.5 ng/uL | 11.0
mCh-RNAi-1 -> 70.2 ng/uL | 7.1
mCh-RNAi-2 -> 75.7 ng/uL | 6.6
mCh-RNAi-3 -> 53.7 ng/uL | 9.3
Total RNA was pooled by treatment so that each sample contributed 1/3rd of its RNA for MinIon sequencing and then speedvacced

#Things to do
differentially expressed gene tables
differentially expressed transcripts tables
htseq-count
deseq2
GO for differentially expressed genes
get bioanalyzer machine and kit specs.


#RANDO
35081 unique transripts in GTF file
34559 unique transcrits in raw ballgown output
FBgn0036091 is a major outlier in bamM male sample, maybe interesting to look at but highly skews gene data, will remove from further analysis

