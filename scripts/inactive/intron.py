#!/usr/bin/env python
import sys, fasta

df = open(sys.argv[1])

num = 0
length = 0
for ident, seq in fasta.FASTAReader(df):
    length += len(seq)
    num += 1
average = length / num
print average