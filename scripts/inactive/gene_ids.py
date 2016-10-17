#!/usr/bin/env python
import sys

try:
    mode = sys.argv[1] #gtf or error
except IndexError:
    mode = ""
stringtie = open("/Users/cmdb/bam_sxl_analysis/data_tables/stringtie.gtf")
gtf = open("/Users/cmdb/bam_sxl_analysis/data_tables/stringtie_updated.gtf")
gene_names = open("/Users/cmdb/annotations/Dmel_gene_names.txt")

ref_gene_dict = {}
for line in gene_names:
    line = line.rstrip("\r\n").split("\t")
    gene_id = line[1]
    gene_name = line[0]
    ref_gene_dict[gene_id] = gene_name

str_gene_dict = {}
for line in stringtie:
    if line.startswith("#"):
        continue
    line = line.rstrip("\r\n").split("\t")
    group = line[8].split(" ")
    if line[1] == "FlyBase":
        continue
    gene_id = group[1]
    if line[2] == "transcript":
        if len(group) == 5:
            gene_ref = gene_id
            if gene_id not in str_gene_dict:
                str_gene_dict[gene_id] = [gene_ref]
            else:
                if gene_ref not in str_gene_dict[gene_id]:
                    str_gene_dict[gene_id].append(gene_ref)
                    #print "gene " + gene_id + " has multiple gene entries " + ", ".join(str_gene_dict[gene_id])
        elif len(group) == 7:
            gene_ref = group[5]
            if gene_id not in str_gene_dict:
                str_gene_dict[gene_id] = [gene_ref]
            else:
                if gene_ref not in str_gene_dict[gene_id]:
                    str_gene_dict[gene_id].append(gene_ref)
                    #print "gene " + gene_id + " has multiple gene entries " + ", ".join(str_gene_dict[gene_id])
    elif line[2] == "exon":
        if len(group) == 7:
            gene_ref = gene_id
            if gene_id not in str_gene_dict:
                str_gene_dict[gene_id] = [gene_ref]
            else:
                if gene_ref not in str_gene_dict[gene_id]:
                    str_gene_dict[gene_id].append(gene_ref)
                    #print "gene " + gene_id + " has multiple gene entries " + ", ".join(str_gene_dict[gene_id])
        if len(group) == 9:
            gene_ref = group[7]
            if gene_id not in str_gene_dict:
                str_gene_dict[gene_id] = [gene_ref]
            else:
                if gene_ref not in str_gene_dict[gene_id]:
                    str_gene_dict[gene_id].append(gene_ref)
                    #print "gene " + gene_id + " has multiple gene entries " + ", ".join(str_gene_dict[gene_id])
str_gene_dict_multi = {}
for key in str_gene_dict:
    if len(str_gene_dict[key]) > 1:
        for i in str_gene_dict[key]:
            if i.startswith("\"MSTRG"):
                continue
            else:
                if key not in str_gene_dict_multi:
                    str_gene_dict_multi[key] = [i]
                else:
                    str_gene_dict_multi[key].append(i)
# for key in str_gene_dict:
#     if len(str_gene_dict[key]) > 1:
#         print str_gene_dict[key]
#         print str_gene_dict_multi[key]

for line in gtf:
    line = line.rstrip("\r\n").split("\t")
    group = line[8].split(" ")
    line = line[:8]
    if group[1].startswith("\"FBgn"):
        print "\t".join(line) + "\t" + " ".join(group)
        continue
    if group[1].startswith("\"MSTR"):
        if group[1] not in str_gene_dict_multi:
            print "\t".join(line) + "\t" + " ".join(group)
            continue
        if len(str_gene_dict_multi[group[1]]) != 1:
            print "\t".join(line) + "\t" + " ".join(group)
            continue
        group[1] = str_gene_dict_multi[group[1]][0]
        gene_ref = group[1]
        gene_ref = gene_ref.replace(";", "").replace("\"", "")
        gene_ref = ref_gene_dict[gene_ref]
        gene_ref = "\"" + gene_ref + ";\""
        group[3] = gene_ref
        print "\t".join(line) + "\t" + " ".join(group)
        