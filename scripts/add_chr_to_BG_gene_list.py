#!/usr/bin/env python
import sys, fasta
import pandas as pd

gene_df = open(sys.argv[1]) #BG gene list object

gtf = open("/Users/cmdb/annotations/r6.12_dmel_noERCC.gtf") #reference gtf file
gtf_dict = {}
for line in gtf:
    line = line.rstrip("\r\n").split("\t")
    group = line[8].split(" ")
    gene = group[1] #gene name
    if gene not in gtf_dict:
        gtf_dict[gene] = [line[0], line[6]]

for i, line in enumerate(gene_df):
    line = line.rstrip( "\r\n").split("\t")
    if i == 0:
        print "chr" + "\t" + "strand" + "\t" + "\t".join(line)
        continue
    gene = "\"" + line[0] + "\";"
    chrom = gtf_dict[gene][0]
    strand = gtf_dict[gene][1]
    line.insert(0, strand)
    line.insert(0, chrom)
    print "\t".join(line)

