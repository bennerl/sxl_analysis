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

def motif_find_exon(target, up_start, down_start, prime):
    for match in re.finditer(motif, target):
            gene = row[13]
            length = match.span()
            motif_match = match.group()
            if row[2] == "+":
                start = up_start + length[0]
                end = up_start + length[1]
            elif row[2] == "-":
                start = down_start - length[1]
                end = down_start - length[0]
            if sys.argv[1] == "bed":
                if prime == "Five_Amb" or prime == "Three_Amb":
                    print "%s%s\t%d\t%d\t%s-%s\t%d\t%s\t%d\t%d\t%s" % ("chr", ident, start, end, gene, prime, 0, row[2], start, end, "0,0,0")
                else:
                    print "%s%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s" % ("chr", ident, start, end, gene, 0, row[2], start, end, "0,0,0")    
            if sys.argv[1] == "stats":
                bamf = float(row[5])
                bamm = float(row[6])                
                mch = float(row[7])
                sxl = float(row[8])
                fc_bam = float(row[9])
                pval_bam = float(row[10])
                fc_rnai = float(row[11])
                pval_rnai = float(row[12])
                print "%s\t%s\t%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s" % (gene, gene, motif_match, fc_rnai, pval_rnai, fc_bam, pval_bam, mch, sxl, bamf, bamm, ident, start, end, row[2], int(row[3]), int(row[4]), prime, row[0] )
        
exon_file = open(sys.argv[2])
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
            up_start = int(row[3]) - 25
            up_end = int(row[3]) + 10
            down_start = int(row[4]) - 10
            down_end = int(row[4]) + 25
            up_seq = seq[up_start:up_end]
            down_seq = seq[down_start:down_end]
            if row[2] == "-":
                up_seq = rev_comp(up_seq)
                down_seq = rev_comp(down_seq)
                motif_find_exon(up_seq, up_start, up_end, "Five")
                motif_find_exon(down_seq, down_start, down_end, "Three")
            elif row[2] == "+":
                motif_find_exon(up_seq, up_start, up_end, "Five")
                motif_find_exon(down_seq, down_start, down_end, "Three")
            else:
                row[2] = "+"
                motif_find_exon(up_seq, up_start, up_end, "Five_Amb")
                motif_find_exon(down_seq, down_start, down_end, "Three_Amb")
                row[2] = "-"
                up_seq = rev_comp(up_seq)
                down_seq = rev_comp(down_seq)
                motif_find_exon(up_seq, up_start, up_end, "Five_Amb")
                motif_find_exon(down_seq, down_start, down_end, "Three_Amb")
    except KeyError:
        continue
