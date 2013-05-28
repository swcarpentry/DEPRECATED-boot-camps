#!/usr/bin/env python
from Bio import SeqIO
import sys

if len(sys.argv) != 3 :
     print "wrong number of args"
     sys.exit()
print sys.argv[1], sys.argv[2]
generator = SeqIO.parse(sys.argv[1], "genbank")
outfile = open(sys.argv[2], "w")

outfile.write(generator.next().format("fasta"))
