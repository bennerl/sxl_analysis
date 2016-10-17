#!/usr/bin/env python
import sys, fasta

if sys.argv[2] == "annotate":
    gtf = open("/Users/cmdb/bam_sxl_analysis/data_tables/r6.12_dmel_noERCC.gtf")
    tran_exon_dict = {}
    gene_tran_dict = {}
    for line in gtf:
        line = line.rstrip("\r\n").split("\t")
        group = line[8].split(" ")
        if line[2] == "gene":
            gene = group[1]
            gene_tran_dict[gene] = []
            continue
        if line[2] == "mRNA" or line[2] == "ncRNA":
            tran = group[5]
            tran_exon_dict[tran] = []
            gene_tran_dict[gene].append(tran)
            continue
        if line[2] == "exon":
            exon = line[0] + "_" + line[3] + "_" + line[4]
            tran_exon_dict[tran].append(exon)
    gtf.close()

    df = open(sys.argv[1])
    for i, line in enumerate(df):
        line = line.rstrip("\r\n").split("\t")
        if i == 0:
            line.append("Annot_Support")
            print "\t".join(line)
            continue
        exon = [line[11], line[15], line[16]]
        gene = "\"" + line[1] + "\";"
        prime = line[17]
        if "_" in gene:
            line.append("No")
            print "\t".join(line)
            continue
        tran_list = gene_tran_dict[gene]
        for tran in tran_list:
            exon_list = tran_exon_dict[tran]
            present = "No"
            length = 0
            for exon_tran in exon_list:
                exon_tran = exon_tran.split("_")
                if length == 0:
                    exon_begin = int(exon_tran[1])
                    exon_end = int(exon_tran[2])
                    length = 1
                else:
                    if int(exon_tran[1]) < exon_begin:
                        exon_begin = int(exon_tran[1])
                    if int(exon_tran[2]) > exon_end:
                        exon_end = int(exon_tran[2])
                if exon == exon_tran:
                    present = "Yes"
                if exon[1] == exon_tran[1] and exon[2] != exon_tran[2] and prime == "Five":
                    present = "Yes"
                if exon[1] != exon_tran[1] and exon[2] == exon_tran[2] and prime == "Three":
                    present = "Yes"
            if int(exon[1]) < exon_begin or int(exon[2]) > exon_end:
                present = "Yes"
            if present == "No":
                line.append("Yes")
                line.append(tran.replace("\"","").replace(";",""))
                print "\t".join(line)
                line.pop()
                line.pop()
                annotate = "Yes"
            else:
                annotate = "No"
        if annotate == "No":
            line.append("No")
            line.append("")
            print "\t".join(line)
            
if sys.argv[2] == "unannotate":
    string = open("/Users/cmdb/bam_sxl_analysis/data_tables/stringtie.gtf")
    string_tran_exon_dict = {}
    string_gene_tran_dict = {}
    string_exon_gene_dict = {}
    for line in string:
        line = line.rstrip("\r\n").split("\t")
        group = line[8].split(" ")
        if line[2] == "transcript" and group[1].startswith("\"MSTR") and group[3].startswith("\"MSTR"):
            gene = group[1]
            tran = group[3]
            if gene not in string_gene_tran_dict:
                string_gene_tran_dict[gene] = [tran]
            else:
                string_gene_tran_dict[gene].append(tran)
            string_tran_exon_dict[tran] = []
            continue
        if line[2] == "exon" and group[1].startswith("\"MSTR") and group[3].startswith("\"MSTR"):
            exon = line[0] + "_" + line[3] + "_" + line[4]
            string_tran_exon_dict[tran].append(exon)
            if exon not in string_exon_gene_dict:
                string_exon_gene_dict[exon] = [gene]
            else:
                if gene not in string_exon_gene_dict[exon]:
                    string_exon_gene_dict[exon].append(gene)
            continue
    string.close()

    df = open(sys.argv[1])
    for i, line in enumerate(df):
        line = line.rstrip("\r\n").split("\t")
        if i == 0:
            line.append("Annot_Support")
            print "\t".join(line)
            continue
        exon = [line[11], line[15], line[16]]
        gene = "\"" + line[1] + "\";"
        prime = line[17]
        if "_" in gene:
            line.append("No")
            print "\t".join(line)
            continue
        try:
            gene_list = string_exon_gene_dict["_".join(exon)]
            for string_gene in gene_list:
                tran_list = string_gene_tran_dict[string_gene]
                for tran in tran_list:
                    exon_list = string_tran_exon_dict[tran]
                    present = "No"
                    length = 0
                    for exon_tran in exon_list:
                        exon_tran = exon_tran.split("_")
                        if length == 0:
                            exon_begin = int(exon_tran[1])
                            exon_end = int(exon_tran[2])
                            length = 1
                        else:
                            if int(exon_tran[1]) < exon_begin:
                                exon_begin = int(exon_tran[1])
                            if int(exon_tran[2]) > exon_end:
                                exon_end = int(exon_tran[2])
                        if exon == exon_tran:
                            present = "Yes"
                        if exon[1] == exon_tran[1] and exon[2] != exon_tran[2] and prime == "Five":
                            present = "Yes"
                        if exon[1] != exon_tran[1] and exon[2] == exon_tran[2] and prime == "Three":
                            present = "Yes"
                    if int(exon[1]) < exon_begin or int(exon[2]) > exon_end:
                        present = "Yes"
                    if present == "No":
                        line.append("Yes")
                        line.append(tran.replace("\"","").replace(";",""))
                        print "\t".join(line)
                        line.pop()
                        line.pop()
                        annotate = "Yes"
                    else:
                        annotate = "No"
                if annotate == "No":
                    line.append("No")
                    line.append("")
                    print "\t".join(line)
        except KeyError:
            line.append("No")
            line.append("")
            print "\t".join(line)
