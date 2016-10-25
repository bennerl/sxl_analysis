#!/usr/bin/env python
import sys, fasta

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

df = open("temp.txt")
seq_file = open("/Users/cmdb/annotations/r6.12_dmel_noERCC.fasta")

chr_dict = {}
for line in df:
    line = line.rstrip("\r\n").split("\t")
    if line[1] not in chr_dict:
        chr_dict[line[1]] = [line]
    else:
        chr_dict[line[1]].append(line)
df.close()

for ident, seq in fasta.FASTAReader(seq_file):
    try:
        for exon in chr_dict[ident]:
            if exon[2] == "-":
                seq_up = seq[(int(exon[3]) - 25) : (int(exon[3]) + 10)]
                seq_down = seq[(int(exon[4]) - 10) : (int(exon[4]) + 25)]
                seq_up = rev_comp(seq_up)
                seq_down = rev_comp(seq_down)
                print ">" + exon[1] + "_" + exon[3] + "_" + exon[4] + "_seq-up"
                print seq_up
                print ">" + exon[1] + "_" + exon[3] + "_" + exon[4] + "_seq-down"
                print seq_down
            else:
                seq_up = seq[(int(exon[3]) - 25) : (int(exon[3]) + 10)]
                seq_down = seq[(int(exon[4]) - 10) : (int(exon[4]) + 25)]
                print ">" + exon[1] + "_" + exon[3] + "_" + exon[4] + "_seq-up"
                print seq_up  
                print ">" + exon[1] + "_" + exon[3] + "_" + exon[4] + "_seq-down"
                print seq_down
    except KeyError:
        continue

