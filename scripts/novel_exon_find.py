#!/usr/bin/env python
import sys

if sys.argv[1] == "gtf":
    stringtie = open("/Users/cmdb/bam_sxl_analysis/data_tables/stringtie.gtf")
    gtf = open("/Users/cmdb/bam_sxl_analysis/data_tables/r6.12_dmel_noERCC.gtf")
    gtf_set = set()
    for line in gtf:
        line = line.rstrip("\r\n").split("\t")
        if line[2] != "exon":
            continue
        exon = line[0] + "_" + line[3] + "_" + line[4]
        if exon not in gtf_set:
            gtf_set.add(exon)

    str_gene_dict = {}
    for line in stringtie:
        line = line.rstrip("\r\n").split("\t")
        group = line[8].split(" ")
        if line[1] == "FlyBase" or line[2] == "transcript" or group[3].startswith("\"FBtr"):
            continue
        exon = line[0] + "_" + line[3] + "_" + line[4]
        if exon not in gtf_set:
            print "\t".join(line[:8]) + "\t" + " ".join(group[:2])

if sys.argv[1] == "merge": 
    gtf = open("/Users/cmdb/bam_sxl_analysis/data_tables/novel_exons.txt")
    gtf_dict = {}
    for line in gtf:
        line = line.rstrip("\r\n").split("\t")
        group = line[8].split(" ")
        exon = line[0] + "_" + line[3] + "_" + line[4] #i checked to make sure each one is unique
        gtf_dict[exon] = group[1]
    df = open(sys.argv[2])
    for line in df:
        line = line.rstrip("\r\n").split("\t")
        exon = line[1] + "_" + line[3] + "_" + line[4]
        if exon in gtf_dict:
            print "\t".join(line) + "\t" + gtf_dict[exon].replace(";", "").replace("\"", "")
