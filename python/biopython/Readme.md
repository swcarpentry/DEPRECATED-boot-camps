
##How to program if you must##

If a bioinformatic task is relatively common, it is likely that someone else has written software to do it already.
Format conversions, paired-read merging, quality trimming, recruitment to references, diploid snp calling, haploid
re-sequencing--these are all problems that can be solved by finding out what software purports to do the job, 
figuring out how to install it, getting the data into the right input formats and getting the results out of 
whatever formats the tools write to.  

Depending on the complexity of the task and the ease-of-use and scope of the existing alternatives, it could be 
easier to adapt an existing software library or it could be easier to write the code to do it yourself.

This section describes programming using the Biopython libraries when you must.
In general you will have your own data, you will need to perform operations on it,
you will get some reference data, perform some comparison operation, and then perform
operations on the results of the comparison in view of your hypotheses.

One note: solving a particular data manipulation problem takes a certain amount of time, time that should include the time to 
fix the mistakes and confirm that your code is really is doing what you think it is.   It is dramatically easier to write a program that
you will use once and then throw away than to write a proram that is useful to you again and again.  Ask: do you really
want to solve a particular (easy) format-conversion problem six times, one for each collaboration?  If you invest effort 
in identifying what you need and what parts you will use against and again, you can forget the details of how you solved
this (boring) problem in the first place and direct your time to more interesting things.

##Biopython##
Be sure to bookmark the Biopython Cookbook and Tutorial
http://biopython.org/DIST/docs/tutorial/Tutorial.html and to cite Biopython in publications that use the tools: Biopython: freely available Python tools for computational molecular biology and bioinformatics. http://www.ncbi.nlm.nih.gov/pubmed/19304878?dopt=Abstract

###Get the data--NCBI's EUTILS###
One large data source is NCBI.  NCBI provides an interface to allow automated download of various data products using HTTP GET requests.

The documentation for this interface, called EFETCH, is here:
http://www.ncbi.nlm.nih.gov/books/NBK25499/

Before using EUTILS, we might want to know what kinds of things it can do:
#####Search engine for PUBMED: #####
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=%22Life+with+6000+genes%22&retmax=100
#####Search engine for SRA#####
 http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=SRX015714
#####Search engine for Genbank#####
 http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=GFAJ1

Ok, you probably get the picture.  These queries return us lists of IDs that don't mean anything to us as humans, but that we can iterate over and retrieve automatically.

Once we've got lists of identifiers, we can retrieve data:
##### Pubmed abstracts:#####
This abstract http://www.ncbi.nlm.nih.gov/pubmed/8849441
can be retrieved in a machine-friendly format with this query:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=8849441&rettype=xml
##### SRA metadata bundles#####
This page of metadata about a dataset from Jeff Gordon's twin study  http://www.ncbi.nlm.nih.gov/sra/SRX015714 
can be retrieved from here:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=16175
##### Genbank records #####
Sequence data deposited by authors (with annotations, when provided by authors or the archive itself) can be retrieved in FASTA or genbank formats using the EUTILS suite.  The following URL should retrieve the human mitochondrial reference sequence (NC_012920) in genbank format:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=251831106&rettype=gb
##### Genome sequences #####
The following query will retrieve the genome of Candidatus Hodgkinia cicadicola (REFseq accession NC_012960.1) in fasta format:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NC_012960.1&rettype=fasta

A table describing the supported formats is here:
http://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.chapter4_table1/?report=objectonly
Note that these queries return data in several different formats, some return XML-formatted data structures while
others return files in gb and fasta formats.  Of course, how you handle the data after retrieving
it is going to depend on the data format.

Biopython has subroutines that take EFETCH's options and parameters and return python objects.
This saves us from having to write code that directly talks to NCBI's EFETCH API, freeing us to spend our time elsewhere.  
We just need to find out how to use these subroutines.

Here is a example using the ```Entrez.efetch```  biopython binding
```python
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

#  Always tell NCBI who you are.  Email address is a mandatory argument.
Entrez.email = "trimble@example.com"
Entrez.tool = "SoftwareCarpentryBootcamp"

# accession = "NC_000913"   # E. coli K12 reference genome accession number

accession = sys.argv[1]     # take the first program argument

downloadstuff(accession)   
```

Exercise:
```ladyslipperITSaccessionnumbers.txt``` contains 94 accession numbers for the ITS ribosomal marker sequences of certain lady slipper orchids.
See if you can modify the recipe above to use FASTA format (instead of genbank format) and then find a way to download all 94 sequences.  Once you have the sequences, you can concatenate them and run your favorite multiple-sequence-alignment program.

#### Ok, just give me the data####
Exercise:
```retrievegenbankfromNCBIbyacession.py``` was intended to accept an accession number on the command line and create a genbank file with the results of the query to NCBI's efetch.  To write the subroutine that gets the data from you will need to study the documentation XXX.
Use it to retrieve NC_001422, the genome of Phi-x174, in genbank format.  We will use this below.

### Iterating through data records ###

SeqIO provides a variety of methods for stepping through data sources one record at a time.  There is variety in the places we
can get the data from, the data types, and the procedures used to access the data.  
Data sources can be web interfaces, filenames, or file handles.  Data types can include amino acid sequences, 
short-read nucleic acid sequences with or without qualities, draft genomes in hundreds of contigs,  or 
complete genomes with gene coordinates, translations, and additional notes about how the genes were identified.  
The access procedures include getting the data from HTTP requests at NCBI, opening a file on disk and loading it all into memory at once, or 
reading a data file one record at a time.  

The SeqIO.parse method returns a *generator*, an python object that will let us see one record at a time.
```
from Bio import SeqIO
generator = SeqIO.parse("NC_000913.gbk", "genbank")
for seqrecord in generator:
    print repr(seqrecord)

```
Aside: The ```Seq``` sequences are by default *immutable*; you will need to use the ```MutableSeq``` data type if you need to actually edit the sequences.

First we describe the data formats, show examples of parsers for each format, and then give some exercises for operating on sequences given a loop over all the sequences in a data source.

Finally, the manual page on SeqIO
```
help(SeqIO)
```
gives recipes for turning sequences into lists 

####Fasta sequence parsing####
The minimal data type for sequence data, this format includes only a text record description and a (possibly long) sequence.  Nucleic acid sequences, partially-ambiguous nucleic acid sequences, and amino acid sequences can all be encoded in this bare-bones format.  

Exercise:  
Modify the existing program ```fastadecorator.py``` to calculate the length and gc content of fasta sequences and to output fasta sequences with `length=XXX gc=XX.X%` appended to the end of the fasta ID.

Exercise:  
Modify the existing program ```fasta2reversecomplement.py``` to output fasta whose sequences have been reverse-complemented.

Exercise:  
Modify the existing program ```fasta2allsixframes.py``` to read in nucleic acid FASTA and output six sequences, representing the translation of each sequence in all six reading frames.

####Fastq sequence parsing####
FASTQ is the de-facto standard data format for the output of modern high-throughput sequencing machines.
In addition to sequence ID and header, this format includes a quality symbol for each base.

Exercise:  
Write a program that removes bases with the lowest FASTQ quality score (encoded as "B") from the end of the read.  That is, for each input sequence, remove all of the base pair from the end that have "B" for the quality score.

####Genbank sequence parsing####

Unlike FASTA and FASTQ, the Genbank format has many optional fields, and the optional fields may be of variable lengths.  
Consequently, Biopython's genbank parser creates ```SeqRecords``` which have variable-length data types inside them.
```
from Bio import SeqIO
generator = SeqIO.parse("NC_001422.gbk", "genbank")
gbrecord = generator.next()  # This grabs the first record.
```

Now let's look it over:

```
type(gbrecord)
    Bio.SeqRecord.SeqRecord
dir(gbrecord) 
   ...
   'annotations',
   'dbxrefs',
   'description',
   'features',
   'format',
   'id',
   'letter_annotations',
   'lower',
   'name',
   'reverse_complement',
   'seq',
   'upper'] 
```

Of these, the attribudes `id` `name` and `description` are top-level descriptions of each sequence (or contig).  These are accessed by `gbrecord.name` `gbrecord.id` and `gb.description`. `grecord.seq` contains the sequence as a `Seq` object. 
`gbrecord.annotations` is a *dict* containing additional details (submitters, date, organism taxonomy) for the sequence as a whole.   
`gbrecord.features` is a *list* of ```SeqFeature``` objects, if provided, that include the genes, protein coding and RNAs, that the creator of the genbank record have decided to include.

Using ```NC_001422.1.gbk``` as an example, let us examine the data structures that we get from parsing it.

```
gbrecord.features
```
The first item in this list is a ```SeqFeature``` object:
```
type(gbrecord.features[0])
   Bio.SeqFeature.SeqFeature
```

This is a list, so we can iterate through it:
```
for i in gbrecord.features:
    print i
```
This shows us that the information is there -- looping through all the features of the gbrecord gives us access to everything except the top-level ```seq``` and top-level ```annotations```, for which there is only one per SeqRecord.  

The ```SeqFeature``` data type has attributes that include  ```id``` ```location``` and ```strand```.  The ```type``` attribute encodes whether the feature represents a gene, tRNA, RNA, or other types, and the available data will usually depend on the value of type.  (RNAs do not have translations; some details are included under "gene" and others under the corresponding "CDS".  A final attribute of ```gbrecord.features``` is ```qualifiers``` -- a *dict* of *list*s of additional annotation fields.

```
for i in gbrecord.features:
    if i.type == "CDS" :
        print i.qualifiers
```
```gbrecord.features[2]``` is the first protein-coding annotation "CDS".  Examining it, we see
```
print gbrecord.features[2].qualifiers
print gbrecord.features[2].qualifiers.keys()
    ['function', 'locus_tag', 'codon_start', 'product', 'transl_table', 'note', 'db_xref', 'translation', 'gene', 'protein_id']

```
We can output a table of all the CDS features, then, by looping over all the features, testing for CDS, and printing the fields we like:

```
for i in gbrecord.features:
    if i.type == "CDS" :
        print i.qualifiers["locus_tag"][0], i.qualifiers["protein_id"][0] , i.qualifiers["gene"][0],  i.qualifiers["translation"][0] 

```
Oops.  This doesn't work.  Not all the CDS records have ```qualifiers["gene"]``` defined, so I get a KeyError when I try to access qualifiers["gene"] for the feature that doesn't have the key "gene".  We need a ```try... except``` conditional to trap the error and fill it in with a default value.
```
for i in gbrecord.features:
    if i.type == "CDS" :
        try: 
            gene = i.qualifiers["gene"][0]
        except KeyError:
            gene = "-"  
        print i.qualifiers["locus_tag"][0], i.qualifiers["protein_id"][0] , gene,  i.qualifiers["translation"][0]

```
Exercise: 
Parse ```NC_001422.1.gbk``` and output a table containing the gene name, the protein ID number, and the translated sequence.  

####XML parsing####
XML is a general-purpose format for structured data; it is a datatype used by some of NCBI's EUTILS and as the most complete machine-readable output format of NCBI BLAST alignments.    

### Similarity searching ###
Similarity searching is perhaps the fundamental operation of computational biology; comparisons between known sequences and novel sequences are the bread-and-butter of sequence interpretation.  
You can run BLAST on your laptop, on your department's server, or via the NCBI web interface.  
BLAST against large databases is an expensive operation, so if your computational plan requires running BLAST a million times, you probably need to re-think your plan.  
For small numbers of sequences and for high-value seuquences (contigs, genomes)  BLAST is extremely popular.


```
    from Bio.Blast import NCBIWWW
    blastresults = NCBIWWW.qblast("blastn", "nr", sequence="CTAAGCACTTGTCTCCTGTTTACTCCCCTGAGCTTGAGGGGTTAACATGAAGGTCATCGATAGCAGGATAATAATACAGTA" )

    result_handle = open("my_blast.xml")

    from Bio.Blast import NCBIXML
    blast_record = NCBIXML.read(result_handle)```

```
    from Bio.Blast.Applications import NcbiblastxCommandline
    help(NcbiblastxCommandline)```

###High-throughput data--getting it###
High-throughput sequencing datasets range in size from a few megabytes to a few hundreds of gigabytes in size.  
Some institutions make raw sequence data available by FTP, but the sequence archive is the largest sequence data warehouse.

The NCBI offers a guide to downloading data here http://www.ncbi.nlm.nih.gov/books/NBK47540/
which includes links to downloading the *SRA toolkit*: http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=std

The sequence read archive maintains its own formats, and its own libary of programs to get data out of the SRA format.  The options for the utilities (and the formats themselves) change from time to time, so if something doesn't work, the first thing the help desk will ask you to do is update your copy of 

PhiX control lane, described at:
http://www.ncbi.nlm.nih.gov/sra/SRX017204
SRR036919
We can download from here 
wget ftp://ftp.ncbi.nih.gov/sra/sra-instant/reads/ByRun/litesra/SRR/SRR036/SRR036919/SRR036919.sra 


