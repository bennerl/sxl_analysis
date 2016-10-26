#Workflow
The workflow that I've completed so far is outlined below.

For all samples, abbreviations are as follows:
BamF -> Bam- females
BamM -> Bam- males
mCh -> Bam- females nanos-gal4 > mCh-RNAi
Sxl -> Bam- females nanos-gal4 > Sxl-RNAi

pos -> means that there was an increase in Sxl sample
neg -> means that there was an increase in mCh sample
FC -> fold change

Reads mapped with Hisat2 ->
Generate guided splicing/gene annotation alignement files (GTF) for each sample with stringtie ->
Merge each samples GTF file with stringtie merge to generate an annotation file with characterized and novel gene/splicesites ->
Generate read alignements and transcript/gene abundance files with stringtie to use with ballgown ->
Determine differential gene expression and exon abundance with ballgown.
(sxl_analysis/scripts/biowulf.sh)

Motif = [A|G][T]{7,}.{0,5}AG.
Motif selection for SXL is somewhat biased, but unfortunately I think is the best I can do with the current methods.  Sakashita & Sakamoto, 1994 did a SELEX-seq experiment for SXL binding to RNA sequences.  They published the motif A(U)(n)(N)(n)AGU.  This is actually not really representative in my opinion though (ie this would not actually call the tra splicesite!), and thus I added the condition that the first base could be A or G (A|G, grabs tra site).  The poly U tract seems pretty legitimate in their data and it seems to be greater than 7, so I put it had to be 7 or more U's to call. THey then claim it has N(n), but realistically their top targets have anywhere from 0 to 3 NT from the poly U to the AG (which is the canonical splicesite).  So I compromised at anywhere from 0 to 5 N nucleotides (some of their next target sequences were up to 5).  THen the AG was totally necessary so I kept that.  Now I searched the end of exons for these sequences, but it probably isnt enriched for this site due to the AG at the end of the motif and thus will be preferential for the 5` splicesite of the exon, but thats purely speculation, havent actually measured that.  There are other motifs/sequences published from crosslinking data.  They are generally much longer and for all published specific motifs, except for one, this motif structure would catch them, so I believe it to be pretty exhaustive for potential SXL motifs at these positional sites.  The one exception is a motif that starts with a C at the beginning of the poly U tract, but I'm not too worried about this one, I can always go back and include it.  

Ballgown analysis was done with FlyBase gtf seperately from novel gtf (ballgown, sxl_data/stringtie_annotated|merge).  The pipeline I used to select for potential targets are visually outlined in sxl_data/powerpoints/pipeline.pptx
Basically, I performed ballgown analysis at the gene and exon level in respect to the FlyBase annotation (script, sxl_analysis/scripts/ballgown_generator_annotated.R; data, sxl_data/data_tables/raw_gene|exon_table_annotated.txt). I used relatively lax requirements for gene and exon selection.  For both, I used a pvalue of < 0.5 (0.5, not 0.05), a gene level FPKM > 1 and exon count level > 10.  I then separated genes into positive or negative fold change (FC) lists for further analysis (script, sxl_analysis/scripts/exon_gene_select.R; data, sxl_data/data_tables/gene|exon_pos|neg_fc_annot.txt). For each set of genes (pos and neg FC), I searched for a modified Sakashita SXL binding site in the 5` and 3` UTR of each gene.  Genes with a motif are then returned and uniquely sorted by FC (raw|uniq_sort, sxl_data/data_tables/Sakashita_gene_pos|neg_annot_stats|.sort_uniq.txt). 

For exon analysis, I extracted -25bps (into intron) and +10bps (into exon) from the 5` splicesite and the reciprocal for the 3` splicesite (+25/-10) relative to the coding strand. I search for a modified Sakashita SXL binding sites within these sequences and return exons with a motif and uniquely sort by FC (raw|uniq_sort, sxl_data/data_tables/Sakashita_exon_pos|neg_annot_stats|.sort_uniq.txt), these are thus potential SXL dependent alternatively spliced exons.  

I then proceed with the rationale that in either a pos or neg exon usage model (pos means increased usage of exon in Sxl-RNAi dataset, theoretically implying that removal of SXL would allow for splicing and use of exon, that would otherwise be skipped in the presence of SXL, opposite rationale for neg), exons would be skipped and this should be supported in the annotations (ie, transcripts with and without this exon).  I therefore take the potential alternatively spliced exons and determine if the FlyBase GTF supports an exon skipping/alternative splicing event for this exon, reporting positive hits (sxl_data/data_tables/Sakashita_exon_pos|neg_annot_support.txt). Exons that did not have support for exon skipping/alternative splicing in the FlyBase annotation were then compared with the novel annotation for novel exon skipping/alternative splicing support and positive hits were reported (sxl_data/data_tables/Sakashita_exon_pos|neg_unannot_support.txt).

I also generated a novel annotation based on the four RNAseq datasets (bamF, bamM, mCh, Sxl) and guided with the FlyBase annotation (sxl_data/data_tables/stringtie.gtf).  I then selected for novel exons in this annotation file (not found in FlYBase Annotation, sxl_data/data_tales/novel_exons.txt) and extracted pos and neg FC exons as described above (sxl_data/data_tables/novel_exon_pos|neg_unannot.txt).  I then repeated the pipeline above selecting for exons with modified Sakashita SXL motifs and were uniquely sorted (sxl_data/data_tables/Sakashita_novel_exon_pos|neg_unnanot_stats|.sort_uniq.txt) and then compared to the novel annotation for exon skipping/alternative splicing support (sxl_data/data_tables/Sakashita_novel_exon_pos|neg_unnanot_supported.txt).

I've added gene names to some of the unprocessed data tables.  These will be named the same as described above except they will contain the name_ extension before them.  So if you open up a table and dont see gene names, look for the same table with the name_ extension in front of the name.

The code for the pipeline is in sxl_analysis/scripts/motif_search_pipeline.sh and is dependent on sxl_analysis/scripts/ add_chr_to_BG_gene_list.py, exon_annotation_compare.py, exon_gene_select.R, exon_unannotated_gene_select.R, fasta.py, novel_exon_find.py, target_finder.py, target_novel_exon.py

I have created a genome browser session for all this data. http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=lbenner&hgS_otherUserSessionName=Sxl
There are three tracks for SXL motifs found with the pipeline described above.  One track for gene (5` and 3` UTR), one track for FlyBase annotated exons, and a third track for novel exons.  The bars represent the SXL motif that was found, green denoting if the motif/gene/exon belong to the pos dataset (increase in exon/gene in Sxl RNAi dataset) and red denoting the neg, opposite direction.  I also supplied the novel guided annotation from stringtie as a track to look at some of the novel exon usage.  Normalized read coverage tracks for each sample are shown in the browser as well. The inverse of the number of reads over 100E6 were multiplied times the raw read depth at each position for library size normalization.  All track files are in the sxl_data/bed/ folder.

#Plots
+++++++++++ I've only generated gene level plots at this time +++++++++++++++

I generated plots of sample FPKM values measuring correlation between individual samples.
MA plots for each genotype/treatment. Specifically look at bamM_FPKM-bamF_FPKM.png and mCh_FPKM-Sxl_FPKM.png
Sxl gene level FPKM and exon count and coverage plot for each sample. sxl_expression.png, sxl_male_exon_count.png, sxl_male_exon_coverage.png
Gene level hierarchical clustering with ward.d and ward.d2 method. clust_gene_ward.d.png and clust_gene_ward.d2.png

Graphing script in sxl_analysis/scripts/graphing.R

#Github
https://github.com/bennerl/sxl_analysis

#MinIon
+++++ this is far from complete, I need to add notes about prep and sequencing run, and need to complete analysis at gene/exon level ++++

Total RNA from above samples were bioanalyzed according to standard protocols.  rRNA bands look great in RNAi samples, a little degraded in bam samples, but not that bad (i've sequenced much worse).  Measured RNA conc. with Qubit HS kit. Sample conc. are as follows:

Sample     -> Conc.      | uL pooled
Sxl-RNAi-1 -> 43.2 ng/uL | 11.6
Sxl-RNAi-2 -> 36.6 ng/uL | 13.7
Sxl-RNAi-3 -> 45.5 ng/uL | 11.0
mCh-RNAi-1 -> 70.2 ng/uL | 7.1
mCh-RNAi-2 -> 75.7 ng/uL | 6.6
mCh-RNAi-3 -> 53.7 ng/uL | 9.3
Total RNA was pooled by treatment so that each sample contributed 1/3rd of its RNA for MinIon sequencing and then speedvacced

have to finish writing up notes on the minion results

#Things to do
exon level plots that correspond to gene level
novel genes that were differentially expressed?
ask caitlin for pdf
ask caitlin for the bioanalyzer details
high end pos and neg FC genes and complete a more exhaustive search for SXl motifs and determine new targets
#RANDO
