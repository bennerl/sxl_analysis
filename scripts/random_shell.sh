################################################################################
rm ~/bam_sxl_analysis/data_tables/stringtie_FBgn.gtf
~/bam_sxl_analysis/scripts/gtf_parse.py > ~/bam_sxl_analysis/data_tables/stringtie_FBgn.gtf

rm ~/bam_sxl_analysis/data_tables/r6.12_dmel_noERCC.browser.gtf
echo "track name=\"FlyBaseAnnotation\" description=\"FlyBase Annotation\" visibility=2 color=205,0,0" > ~/bam_sxl_analysis/data_tables/r6.12_dmel_noERCC.browser.gtf
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9" "$10" "$11" "$12" "$13" "$14" "$15" "$16}' ~/bam_sxl_analysis/data_tables/r6.12_dmel_noERCC.gtf | grep -v chrmitochondrion_genome | grep -v chrUnmapped_Scaffold | grep -v chr2110000 | grep -v chrrDNA | grep -v FBtr0307760 | grep -v FBtr0307759 | grep -v FBtr0084085 | grep -v FBtr0084084 | grep -v FBtr0084083 | grep -v FBtr0084082 | grep -v FBtr0084081 | grep -v FBtr0084080 | grep -v FBtr0084079 >> ~/bam_sxl_analysis/data_tables/r6.12_dmel_noERCC.browser.gtf

rm ~/bam_sxl_analysis/data_tables/stringtie.browser.gtf
echo "track name=\"StringtieAnnotation\" description=\"Stringtie Annotation\" visibility=2 color=0,0,0" > ~/bam_sxl_analysis/data_tables/stringtie.browser.gtf
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9" "$10" "$11" "$12" "$13" "$14" "$15" "$16}' ~/bam_sxl_analysis/data_tables/stringtie_FBgn.gtf | grep -v chrmitochondrion_genome | grep -v chrUnmapped_Scaffold | grep -v chr2110000 | grep -v chrrDNA | grep -v FBtr0307760 | grep -v FBtr0307759 | grep -v FBtr0084085 | grep -v FBtr0084084 | grep -v FBtr0084083 | grep -v FBtr0084082 | grep -v FBtr0084081 | grep -v FBtr0084080 | grep -v FBtr0084079 | grep -v chr3Cen_mapped_Scaffold | grep -v chrX3X4_mapped_Scaffold | grep -v chrXY_mapped_Scaffold | grep -v chrY_mapped_Scaffold >> ~/bam_sxl_analysis/data_tables/stringtie.browser.gtf
################################################################################
