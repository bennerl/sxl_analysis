#!/usr/bin/env python
import sys


stringtie = open("/Users/cmdb/bam_sxl_analysis/data_tables/stringtie.gtf")
gene_names = open("/Users/cmdb/annotations/Dmel_gene_names.txt")
gene_id = ""

ref_gene_dict = {}
for line in gene_names:
    gene_names = open("/Users/cmdb/annotations/Dmel_gene_names.txt")
    line = line.rstrip("\r\n").split("\t")
    gene_id = line[1]
    gene_name = line[0]
    ref_gene_dict[gene_id] = gene_name
gene_names.close()

for line in stringtie:
    if line.startswith("#"):
        continue
    line = line.rstrip("\r\n").split("\t")
    group = line[8].split(" ")
    line = line[:8]
    group.pop()
    if line[1] == "FlyBase":
        gene_ref = group[1]
        gene_ref_strip = gene_ref.replace(";", "").replace("\"", "")
        gene_name = ref_gene_dict[gene_ref_strip]
        gene_name = "\"" + gene_name + ";\""
        group.pop()
        group.insert(2, "gene_symbol")
        group.insert(3, gene_name)
        print "\t".join(line) + "\t" + " ".join(group)
        continue        
    temp_gene = group[1]
    if temp_gene != gene_id:
        gene_id = group[1]
        try:
            gene_ref = group[5]
            gene_ref_strip = gene_ref.replace(";", "").replace("\"", "")
            gene_name = ref_gene_dict[gene_ref_strip]
            gene_name = "\"" + gene_name + ";\""
        except IndexError:
            gene_ref = gene_id
            gene_name = gene_id
        if line[2] == "transcript":
            if len(group) == 6:
                group[1] = gene_ref
                group.pop()
                group.pop()
                group.insert(2, "gene_symbol")
                group.insert(3, gene_name)
                print "\t".join(line) + "\t" + " ".join(group)
            if len(group) == 4:
                group[1] = gene_ref
                group.insert(2, "gene_symbol")
                group.insert(3, gene_name)
                print "\t".join(line) + "\t" + " ".join(group)
        elif line[2] == "exon":
            if len(group) == 8:
                group[1] = gene_ref
                group.pop()
                group.pop()
                group.insert(2, "gene_symbol")
                group.insert(3, gene_name)
                print "\t".join(line) + "\t" + " ".join(group)
            if len(group) == 6:
                group[1] = gene_ref
                group.insert(2, "gene_symbol")
                group.insert(3, gene_name)
                print "\t".join(line) + "\t" + " ".join(group)
    else:
        if line[2] == "transcript":
            if len(group) == 6:
                group[1] = gene_ref
                group.pop()
                group.pop()
                group.insert(2, "gene_symbol")
                group.insert(3, gene_name)
                print "\t".join(line) + "\t" + " ".join(group)
            if len(group) == 4:
                group[1] = gene_ref
                group.insert(2, "gene_symbol")
                group.insert(3, gene_name)
                print "\t".join(line) + "\t" + " ".join(group)
        elif line[2] == "exon":
            if len(group) == 8:
                group[1] = gene_ref
                group.pop()
                group.pop()
                group.insert(2, "gene_symbol")
                group.insert(3, gene_name)
                print "\t".join(line) + "\t" + " ".join(group)
            if len(group) == 6:
                group[1] = gene_ref
                group.insert(2, "gene_symbol")
                group.insert(3, gene_name)
                print "\t".join(line) + "\t" + " ".join(group)
