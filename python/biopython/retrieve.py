#!/usr/bin/env python
import os,sys
from Bio import Entrez
from Bio import SeqIO

def downloadstuff(accessionno):
    filename = "%s.gbk" % accessionno       
    print "Trying efectch on %s, writing to %s" % ( accessionno, filename )
    if not os.path.isfile(filename):  
        net_handle = Entrez.efetch(db="nucleotide",id=accessionno,rettype="gb", retmode="text") 
        out_handle = open(filename, "w")
        out_handle.write(net_handle.read() ) 
        out_handle.close()
        net_handle.close()
    else:
        print "skipping, %s already exists!" % filename

Entrez.email = "trimble@example.com"
Entrez.tool = "SoftwareCarpentryBootcamp"

# accession = "NC_000913"   # E. coli K12 reference genome accession number

accession = sys.argv[1]     # take the first program argument

downloadstuff(accession)    
