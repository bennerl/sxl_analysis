#!/usr/bin/env python
import sys, fasta, re

motif = '[A|G][T]{7,}.{0,5}AG.'
seq_file = open("/Users/cmdb/annotations/r6.12_dmel_noERCC.fasta")
table_type = sys.argv[1]

def rev_comp(sequence):
    revcomp = []
    for i in reversed(sequence):
        if i == "A":
            revcomp.append("T")
        elif i == "T":
            revcomp.append("A")
        elif i == "G":
            revcomp.append("C")
        elif i == "C":
            revcomp.append("G")
    return "".join(revcomp)

def gene_names_find(gene):
    if "_" in gene:
        gene = gene.split("_")
        gene_id = []
        for i in gene:
            gene_id.append(gene_names_dic[i])
        return "_".join(gene_id)
    else:
        return gene_names_dic[gene]

def motif_find_exon(target, up_start, down_start, prime):
    for match in re.finditer(motif, target):
            gene = e2g_dict[e_id]
            if len(gene) > 1:
                gene = "_".join(gene)
            else:
                gene = gene[0]
            length = match.span()
            motif_match = match.group()
            if strand == "+":
                start = up_start + length[0]
                end = up_start + length[1]
            else:
                start = down_start - length[1]
                end = down_start - length[0]
            exon_start = row[3]
            exon_end = row[4]
            if sys.argv[2] == "bed":
                print "%s%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s" % ("chr", ident, start, end, gene, 0, strand, start, end, "0,0,0")
            if sys.argv[2] == "stats":
                bamf = float(row[5])
                bamm = float(row[6])                
                mch = float(row[7])
                sxl = float(row[8])
                fc_rnai = float(row[11])
                pval_rnai = float(row[12])
                fc_bam = float(row[9])
                pval_bam = float(row[10])
                gene_id = gene_names_find(gene)
                print "%s\t%s\t%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s" % (gene_id, gene, motif_match, fc_rnai, pval_rnai, fc_bam, pval_rnai, mch, sxl, bamf, bamm, ident, start, end, strand, exon_start, exon_end, prime,row[0] )

gene_names = open("/Users/cmdb/annotations/Dmel_gene_names.txt")
gene_names_dic = {}
for gene_names_line in gene_names:
    gene_names_line = gene_names_line.rstrip("\r\n").split("\t")
    if gene_names_line[1] not in gene_names_dic:
        gene_names_dic[gene_names_line[1]] = gene_names_line[0]

if table_type == "exon":
    
    e2g = open("/Users/cmdb/sxl_data/data_tables/e2g_annotated.txt")
    e2g_dict = {}
    for e2g_line in e2g:
        e2g_line = e2g_line.rstrip("\r\n").split("\t")
        if e2g_line[0] not in e2g_dict:
            e2g_dict[e2g_line[0]] = [e2g_line[1]]
        else:
            e2g_dict[e2g_line[0]].append(e2g_line[1])
            
    exon_file = open(sys.argv[3])
    exon_file_dict = {}
    for exon_file_line in exon_file:
        exon_file_line = exon_file_line.rstrip("\r\n").split("\t")
        if exon_file_line[1] not in exon_file_dict:
            exon_file_dict[exon_file_line[1]] = [exon_file_line]
        else:
            exon_file_dict[exon_file_line[1]].append(exon_file_line)
    for ident, seq in fasta.FASTAReader(seq_file):
        try:
            for row in exon_file_dict[ident]:
                e_id = row[0]
                strand = row[2]
                up_start = int(row[3]) - 25
                up_end = int(row[3]) + 10
                down_start = int(row[4]) - 10
                down_end = int(row[4]) + 25
                up_seq = seq[up_start:up_end]
                down_seq = seq[down_start:down_end]
                if strand == "-":
                    up_seq = rev_comp(up_seq)
                    down_seq = rev_comp(down_seq)
                    motif_find_exon(up_seq, up_start, up_end, "Five")
                    motif_find_exon(down_seq, down_start, down_end, "Three")
                else:
                    motif_find_exon(up_seq, up_start, up_end, "Five")
                    motif_find_exon(down_seq, down_start, down_end, "Three")
        except KeyError:
            continue

def motif_find_gene(ident, gene, target, location, strand, attribute):
    for match in re.finditer(motif, target):
        length = match.span()
        if strand == "+":
            start = location[0] + length[0]
            end = location[0] + length[1]
        else:
            start = location[1] - length[1]
            end = location[1] - length[0] 
        if sys.argv[2] == "bed":
            print "%s%s\t%d\t%d\t%s-%s\t%d\t%s\t%d\t%d\t%s" % ("chr", ident, start, end, gene, attribute, 0, strand, start, end, "0,0,0")
        if sys.argv[2] == "stats":
            gene_id = gene_names_find(gene)
            motif_match = match.group()
            bamf = float(row_gene_file[3])
            bamm = float(row_gene_file[4])
            mch = float(row_gene_file[5])
            sxl = float(row_gene_file[6])
            fc_bam = float(row_gene_file[7])
            pval_bam = float(row_gene_file[8])
            fc_rnai = float(row_gene_file[9])
            pval_rnai = float(row_gene_file[10])
            print "%s\t%s\t%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%s\t%d\t%d\t%s\t%s\t%s" % (gene_id, gene, motif_match, fc_rnai, pval_rnai, fc_bam, pval_bam, mch, sxl, bamf, bamm, ident, start, end, strand, attribute, transcript.replace(";", "").replace("\"", "") )

if table_type == "gene":
    
    gene_file = open(sys.argv[3])
    gene_file_dic = {}
    for i, gene_file_line in enumerate(gene_file):
        if i == 0:
            continue
        gene_file_line = gene_file_line.rstrip("\r\n").split("\t")
        key = gene_file_line[0]
        if gene_file_line[0] not in gene_file_dic:
            gene_file_dic[gene_file_line[0]] = [gene_file_line]
        else:
            gene_file_dic[gene_file_line[0]].append(gene_file_line)

#dict of chr:[line]
    gtf_file = open("/Users/cmdb/sxl_data/data_tables/r6.12_dmel_noERCC.gtf")
    UTR_dict = {}
    for gtf_file_line in gtf_file:
        gtf_file_line = gtf_file_line.rstrip("\r\n").split("\t")
        if gtf_file_line[2] != "5UTR" and gtf_file_line[2] != "3UTR":
            continue
        gtf_file_gene = gtf_file_line[8].split(" ")
        gtf_file_gene = gtf_file_gene[1].replace("\"", "").replace(";", "")
        if gtf_file_gene not in UTR_dict:
            UTR_dict[gtf_file_gene] = [gtf_file_line]
        else:
            UTR_dict[gtf_file_gene].append(gtf_file_line)
#dict of gene:[UTR-line]
    for ident, seq in fasta.FASTAReader(seq_file):
        try:
            for row_gene_file in gene_file_dic[ident]:
                gene = row_gene_file[2]
                strand = row_gene_file[1]
                try:
                    for row_UTR in UTR_dict[gene]:
                        transcript = row_UTR[8].split(" ")
                        transcript = transcript[5]
                        location = [int(row_UTR[3]), int(row_UTR[4])]
                        target = seq[location[0]:location[1]]
                        attribute = row_UTR[2]
                        if strand == "-":
                            target = rev_comp(target)
                        motif_find_gene(ident, gene, target, location, strand, attribute)
                except KeyError:
                    continue
        except KeyError:
            continue