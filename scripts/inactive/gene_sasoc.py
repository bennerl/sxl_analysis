#!/usr/bin/env python
import sys

stringtie = open("/Users/cmdb/bam_sxl_analysis/data_tables/stringtie.gtf")
# gene_names = open("/Users/cmdb/annotations/Dmel_gene_names.txt")
#
# ref_gene_dict = {}
# for line in gene_names:
#     line = line.rstrip("\r\n").split("\t")
#     gene_id = line[1]
#     gene_name = line[0]
#     ref_gene_dict[gene_id] = gene_name


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
        print "\t".join(line[:8]) + "\t" + "".join(group[:2])
    
    
    
    
#
#
#
#     if line[2] == "transcript":
#         continue
#
#     if len(group) == 6:
#         if group[3].startswith("\"FBtr"):
#             continue
#         tran_ref = group[3]
#         if gene_id not in str_gene_dict:
#             str_gene_dict[gene_id] = [tran_ref]
#         else:
#             if tran_ref not in str_gene_dict[gene_id]:
#                 str_gene_dict[gene_id].append(tran_ref)
#     else:
#         if group[3].startswith("\"FBtr"):
#             continue
#         tran_ref = group[3]
#         if gene_id not in str_gene_dict:
#             str_gene_dict[gene_id] = [tran_ref]
#         else:
#             if tran_ref not in str_gene_dict[gene_id]:
#                 str_gene_dict[gene_id].append(tran_ref)
# for key in str_gene_dict:
#     print key + "\t" + ",".join(str_gene_dict[key])
#
# stringtie.close()

# MSTRG_dict = {}
# for key in str_gene_dict:
#     if len(str_gene_dict[key]) == 1:
#         if str_gene_dict[key][0].startswith("\"MSTRG"):
#             continue
#         MSTRG_dict[key] = str_gene_dict[key][0]
#     if len(str_gene_dict[key]) == 2:
#         if str_gene_dict[key][0].startswith("\"FBgn") and str_gene_dict[key][1].startswith("\"FBgn"):
#             continue
#         for i in str_gene_dict[key]:
#             if i.startswith("\"MSTRG"):
#                 continue
#             MSTRG_dict[key] = i
#
# stringtie = open("/Users/cmdb/bam_sxl_analysis/data_tables/stringtie.gtf")
#
# for line in stringtie:
#     line = line.rstrip("\r\n").split("\t")
#     if line[1] == "FlyBase":
#         print "\t".join(line)
#         continue
#     group = line[8].split(" ")
#     line = line[:8]
#     group.pop()
#     gene = group[1]
#     if gene in MSTRG_dict:
#         group[1] = MSTRG_dict[gene]
#     print "\t".join(line) + "\t" + " ".join(group)

