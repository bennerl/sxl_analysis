#!/usr/bin/env python
import sys

df = open(sys.argv[1])

for line in df:
    line = line.rstrip("\r\n").split(" ")
    if line.startswith("#"):
        continue
    if line.startswith("a"):
        continue
    if line.startswith("s") and "Basecall_2D" not in line:
        print line[1]
        print line[2]