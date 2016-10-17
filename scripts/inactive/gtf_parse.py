#!/usr/bin/env python
import sys

stringtie = open("/Users/cmdb/bam_sxl_analysis/data_tables/stringtie.gtf")
gene_names = open("/Users/cmdb/annotations/Dmel_gene_names.txt")

ref_gene_dict = {}
for line in gene_names:
    line = line.rstrip("\r\n").split("\t")
    gene_id = line[1]
    gene_name = line[0]
    ref_gene_dict[gene_id] = gene_name

str_gene_dict = {}
for line in stringtie:
    line = line.rstrip("\r\n").split("\t")
    group = line[8].split(" ")
    group.pop()
    if line[1] == "FlyBase":
        continue
    gene_id = group[1]
    if line[2] != "transcript":
        continue
    if len(group) == 6:
        gene_ref = group[5]
        if gene_id not in str_gene_dict:
            str_gene_dict[gene_id] = [gene_ref]
        else:
            if gene_ref in str_gene_dict[gene_id]:
                continue
            else:
                str_gene_dict[gene_id].append(gene_ref)
    else:
        gene_ref = gene_id
        if gene_id not in str_gene_dict:
            str_gene_dict[gene_id] = [gene_ref]
        else:
            if gene_ref in str_gene_dict[gene_id]:
                continue
            else:
                str_gene_dict[gene_id].append(gene_ref)
                
stringtie.close()

MSTRG_dict = {}
for key in str_gene_dict:
    if len(str_gene_dict[key]) == 1:
        if str_gene_dict[key][0].startswith("\"MSTRG"):
            continue
        MSTRG_dict[key] = str_gene_dict[key][0]
    if len(str_gene_dict[key]) == 2:
        if str_gene_dict[key][0].startswith("\"FBgn") and str_gene_dict[key][1].startswith("\"FBgn"):
            continue
        for i in str_gene_dict[key]:
            if i.startswith("\"MSTRG"):
                continue
            MSTRG_dict[key] = i

stringtie = open("/Users/cmdb/bam_sxl_analysis/data_tables/stringtie.gtf")

for line in stringtie:
    line = line.rstrip("\r\n").split("\t")
    if line[1] == "FlyBase":
        print "\t".join(line)
        continue
    group = line[8].split(" ")
    line = line[:8]
    group.pop()
    gene = group[1]
    if gene in MSTRG_dict:
        group[1] = MSTRG_dict[gene]
    print "\t".join(line) + "\t" + " ".join(group)


