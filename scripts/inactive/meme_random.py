#!/usr/bin/env python
import sys, fasta, random

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

num_line = int(sys.argv[1]) - 1 #number of lines in file
num_exon = int(sys.argv[2]) #number of exons to choose
num_set = set()
for i in range(0, num_exon):
    x = 0
    while x == 0:
        temp = random.randint(0,num_line)
        if temp not in num_set:
            num_set.add(temp)
            x = 1
        else:
            x = 0

df = open(sys.argv[3])
seq_file = open("/Users/cmdb/annotations/r6.12_dmel_noERCC.fasta")

chr_dict = {}
for i, line in enumerate(df):
    if i not in num_set:
        continue
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
                seq_search = seq[(int(exon[4]) - 10) : (int(exon[4]) + 25)]
                seq_search = rev_comp(seq_search)
                print ">" + exon[1] + "_" + exon[3] + "_" + exon[4] + "_seq-down"
                print seq_search
            else:
                seq_search = seq[(int(exon[3]) - 25) : (int(exon[3]) + 10)]
                print ">" + exon[1] + "_" + exon[3] + "_" + exon[4] + "_seq-up"
                print seq_search
    except KeyError:
        continue

