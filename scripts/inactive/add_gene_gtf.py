#!/usr/bin/env python
import sys
import pandas as pd

gtf = open("/Users/cmdb/bam_sxl_analysis/data_tables/stringtie_updated_2.gtf")
gene = ""

def gene_span(gene):
    num = 1
    for i in dic[gene]:
        if num == 1:
            start = int(i[3])
            end = int(i[4])
            num += 1
        else:
            if int(i[3]) < start:
                start = int(i[3])
            if int(i[4]) > end:
                end = int(i[4])
    return [start, end]

dic = {}
for line in gtf:
    line = line.rstrip("\r\n").split("\t")
    gene = line[8].split(" ")[1]
    if gene not in dic:
        dic[gene] = [line]
    else:
        dic[gene].append(line)
gtf.close()

gtf = open("/Users/cmdb/bam_sxl_analysis/data_tables/stringtie_updated_2.gtf")
for line in gtf:
    line = line.rstrip("\r\n").split("\t")
    group = line[8].split(" ")
    if line[2] == "transcript":
        temp_gene = group[1]
        if temp_gene != gene:
            gene = group[1]
            location = gene_span(gene)
            print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (line[0], line[1], "gene", location[0], location[1], ".", line[6], line[7], line[8])
            print "\t".join(line)
        else:
            print "\t".join(line)
    else:
        print "\t".join(line)


    



