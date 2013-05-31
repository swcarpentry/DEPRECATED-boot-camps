#!/usr/bin/env python
'''exercis-se_allsixframes.py   Write the code to translate sequences 
in all six frames.  ''' 

from Bio import SeqIO

def allsixframes(record):
    ''' This function takes a SeqRecord object and returns a list of 
    six strings, with the translation of the sequence in the seqrecord 
    in all six frames.'''
#   your code goes here
    frame = []
    frame.append("NNN")  # translation
    frame.append("NNN")  # translation + 1
    frame.append("NNN")  # translation + 2
    frame.append("NNN")  # reverse complement translation 
    frame.append("NNN")  # reverse complement translation + 1 
    frame.append("NNN")  # you should see the pattern by now
    return frame

#   Open a fastq file, goes through it record-by-record, and output
#   the sequence id, the sequence, and the translations 

generator = SeqIO.parse("data/test-sequences.fasta", "fasta")
for seqrecord in generator:
    sixframes = allsixframes(seqrecord)
    print ">%s\n%s" % (seqrecord.id, seqrecord.seq)
    for i in range(6):
        print i, sixframes[i] 
