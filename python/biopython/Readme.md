
##How to program if you must##

If a bioinformatic task is relatively common, it is likely that **someone else has written software to do it already**.

Who doesn't want to push the boring, easy tasks to machines?

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
in identifying what you need and what parts you will use again and again, you can forget the details of how you solved
this (boring) problem in the first place and direct your time to more interesting things.

##Python is useful##
One day, a colleague of mine showed me that the MG-RAST website had an interface that would deliver a bundle of data about a dataset in response to an HTTP request.  
Specifically, the request 
http://api.metagenomics.anl.gov/metagenome_statistics/mgm4440613.3?verbosity=full
Has tables of numbers representing the GC content, the length distribution, and high-level summaries of the taxonomic annotations of an NGS dataset.
The script ```metagenome_statistics-example.py``` contains example code that retrieves the metagenome summary data structure from the website.  This data is in a JSON object.  There is a python module to parse the JSON into a python dict of dict, and and plots a sequence-length distribution cotained in the response.  

##What to expect##
We have to **get the data**, **get the data out of its container**, and **do something with the data**.  

##Biopython##
Be sure to bookmark the Biopython Cookbook and Tutorial
http://biopython.org/DIST/docs/tutorial/Tutorial.html and to cite Biopython in publications that use the tools: *Biopython: freely available Python tools for computational molecular biology and bioinformatics.* Bioinformatics. 2009 Jun 1;25(11):1422-3. doi: 10.1093/bioinformatics/btp163. http://www.ncbi.nlm.nih.gov/pubmed/19304878?dopt=Abstract

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

Here is a example using the ```Entrez.efetch```  biopython methods
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

Entrez.email = "trimble@example.com"        #  Always tell NCBI who you are.  
Entrez.tool = "SoftwareCarpentryBootcamp"

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

The SeqIO.parse method returns a *generator*, an python object that will let us see one record at a time.  The records that it returns are ```SeqRecord``` classes, a data type that can accomodate fasta, fastq, and genbank formatted-data, though with different levels of detail.

We can open ```data/tiny-fasta.fasta``` and report the ```type``` of the first record like this:

```python
from Bio import SeqIO
firstrecord = SeqIO.parse("data/tiny-fasta.fasta", "fasta").next()
print type(firstrecord)
```
Note: In an interactive session, we can omit the ```from Bio import SeqIO```.

And we can loop through all the records, one at a time, like this:
```python
generator = SeqIO.parse("data/tiny-fasta.fasta", "fasta")
for seqrecord in generator:
    print repr(seqrecord)
```

Note: ```SeqIO``` takes the format as a mandatory second parameter.  ```fasta``` ```fastq```  ```genbank``` and ```embl``` are among the more important supported formats.

Aside: The ```Seq``` sequences are by default *immutable*; you won't be able to change them unless you use the ```MutableSeq``` data type.

First we describe the data formats, show examples of parsers for each format, and then give some exercises for operating on sequences given a loop over all the sequences in a data source.

Finally, you can get reminders of the recipes for accessing data and specificaitons for the parameters from the manual page on ```SeqIO```:
```python
help(SeqIO)
```

####Fasta sequence parsing####
The minimal data type for sequence data, this format includes only a text record description and a (possibly long) sequence.  Nucleic acid sequences, partially-ambiguous nucleic acid sequences, and amino acid sequences can all be encoded in this bare-bones format.  

```SeqRecord``` data types have the attributes
+ ```.name```  which is the **fasta id** -- all the text before the first whitespace on the header line
+ ```.description``` the entire header line, including the fasta id and anything after the first whitespace
+ ```.seq```  the sequence, as a ```SeqIO.Seq object```

```python
from Bio import SeqIO
for seqrecord in SeqIO.parse("tiny-fasta.fasta", "fasta"):
     print seqrecord.name, seqrecord.seq
```

Exercise:  
Modify the existing program ```fasta2reversecomplement.py``` to output fasta whose sequences have been reverse-complemented.

Exercise:  
Modify the existing program ```fasta2allsixframes.py``` to read in nucleic acid FASTA and output six sequences, representing the translation of each sequence in all six reading frames.  Hint:  Slicing works on ```Seq``` objects like it does on strings, so if ```seq``` is one of the ```seqrecord.seq``` objects, 
```seq[1:]``` and ```seq[2:]``` the sequence with the first character chopped off, and the sequence with the first two characters chopped off, respectively.  
Hint: ```Seq``` objects also have  ```seq.reverse_complement()``` and ```seq.translate()``` methods that return ```Seq``` objects with the reverse complement and the code-11 translation.

If you want to manipulate (say, output) the ```Seq``` objects yourself, ```str(seqrecord.seq)``` will return a string representation of the sequence.  

####Fastq sequence parsing####
FASTQ is the de-facto standard data format for the output of modern high-throughput sequencing machines.
In addition to sequence ID and header, this format includes a quality symbol for each base.

One way to parse fastq is using exactly the same ```SeqIO.parse()``` method, just with ```fastq``` instead of ```fasta``` as the format parameter.

Another approach is to use ```FastqGeneralIterator```.  Unlike ```SeqIO.parse()``` it takes file handles (not file names) and has no format parameter (it only works for fastq).  It returns tuples with the sequence description, the sequence string, and the quality string without additional methods to interpret and format the results.   The following code snippet opens the file tiny-fastq.fastq and writes truncated versions of the data to standard out.  (Note that if any of the input sequences are less than 30 base pairs in length, this code breaks.) 

```python
from Bio.SeqIO.QualityIO import FastqGeneralIterator
in_handle = open("tiny-fastq.fastq")
iterator = FastqGeneralIterator(in_handle)
for triplet in iterator:
     (description, sequence, quality) = triplet 
     print "@%s\n%s\n+\n%s" % ( description, sequence[0:30], quality[0:30] )
```

Exercise:  
Write a program that removes bases with the lowest FASTQ quality score (encoded as "B") from the end of the read.  That is, for each input sequence, remove all of the base pair from the end that have "B" for the quality score.

####Genbank sequence parsing####

Unlike FASTA and FASTQ, the Genbank format has many optional fields, and the optional fields may be of variable lengths.  
Consequently, Biopython's genbank parser creates ```SeqRecords``` which have variable-length data types inside them.

Using ```NC_001422.1.gbk``` as an example, let us examine the data structures that we get from parsing it.

```python
from Bio import SeqIO
generator = SeqIO.parse("NC_001422.gbk", "genbank")
gbrecord = generator.next()  # This grabs the first record.
```

Now let's look it over:

```python
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

Of these, the attribudes `id` `name` and `description` are top-level attributes which describe each sequence (contig).  These are accessed by `gbrecord.name` `gbrecord.id` and `gbrecord.description`. `grecord.seq` contains the sequence as a `Seq` object. 
`gbrecord.annotations` is a *dict* containing additional details (submitters, date, organism taxonomy) for the sequence as a whole.   
`gbrecord.features` is a *list* of ```SeqFeature``` objects, if provided, that include the genes, protein coding and RNAs, that the creator of the genbank record have decided to include.

```python
gbrecord.features
```
The first item in this list is a ```SeqFeature``` object:
```
type(gbrecord.features[0])
   Bio.SeqFeature.SeqFeature
```

This is a list, so we can iterate through it:
```python
for i in gbrecord.features:
    print i
```
This shows us that the information is there -- looping through all the features of the gbrecord gives us access to everything except the top-level ```seq``` and top-level ```annotations```, for which there is only one per SeqRecord.  

The ```SeqFeature``` data type has attributes that include  ```id``` ```location``` and ```strand```.  The ```type``` attribute encodes whether the feature represents a gene, tRNA, RNA, or other types, and the available data will usually depend on the value of type.  (RNAs do not have translations; some details are included under "gene" and others under the corresponding "CDS".  A final attribute of ```gbrecord.features``` is ```qualifiers``` -- a *dict* of *list*s of additional annotation fields.

```python
for i in gbrecord.features:
    if i.type == "CDS" :
        print i.qualifiers
```
```gbrecord.features[2]``` is the first protein-coding annotation "CDS".  Examining it, we see
```python
print gbrecord.features[2].qualifiers
print gbrecord.features[2].qualifiers.keys()
    ['function', 'locus_tag', 'codon_start', 'product', 'transl_table', 'note', 'db_xref', 'translation', 'gene', 'protein_id']

```
We can output a table of all the CDS features, then, by looping over all the features, testing for CDS, and printing the fields we like:

```python
for i in gbrecord.features:
    if i.type == "CDS" :
        print i.qualifiers["locus_tag"][0], i.qualifiers["protein_id"][0] , i.qualifiers["gene"][0],  i.qualifiers["translation"][0] 

```
Oops.  This doesn't work.  Not all the CDS records have ```qualifiers["gene"]``` defined, so I get a KeyError when I try to access qualifiers["gene"] for the feature that doesn't have the key "gene".  We need a ```try... except``` conditional to trap the error and fill it in with a default value.
```python
for i in gbrecord.features:
    if i.type == "CDS" :
        try: 
            gene = i.qualifiers["gene"][0]
        except KeyError:
            gene = "-"  
        print i.qualifiers["locus_tag"][0], i.qualifiers["protein_id"][0] , gene,  i.qualifiers["translation"][0]

```
Exercise: 
Parse ```NC_001422.1.gbk``` and generate an amino acid fasta file containing the translations of all the coding sequences.  Hint: What field do you want to use for the FASTA ID?  Do you want to put anything else in the fasta description line? 

Parse the ```NC_001422.1.faa``` and generate an amino acid fasta file containing the translations of all the coding sequences.  Hint: What field do you want to use for the FASTA ID?  Do you want to put anything else in the fasta description line? 
 

####XML parsing####
XML is a general-purpose format for structured data; it is a datatype used by some of NCBI's EUTILS and as the most complete machine-readable output format of NCBI BLAST alignments.    

### Similarity searching ###
Similarity searching is perhaps the fundamental operation of computational biology; comparisons between known sequences and novel sequences are the bread-and-butter of sequence interpretation.  
You can run BLAST on your laptop, on your department's server, or via the NCBI web interface.  
BLAST against large databases is an expensive operation, so if your computational plan requires running BLAST a million times, you probably need to re-think your plan.  
For small numbers of sequences and for high-value sequences (contigs, genomes)  BLAST is extremely popular.


```python
    from Bio.Blast import NCBIWWW
    mysterysequence = "GCACTTGTCTCCTGTTTACTCCCCTGAGCTTGAGGGGTTAACATGAAGGTCATCGATAGCAGGATAATAATACAGTA"
    blastresults = NCBIWWW.qblast("blastn", "nr", sequence=mysterysequence ) 

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

The sequence read archive maintains its own formats, and its own libary of programs to get data out of the SRA format.  The options for the utilities (and the formats themselves) change from time to time, so if something doesn't work, the first thing the help desk will ask you to do is update your copy of the sra toolkit.
wget ftp://ftp.ncbi.nih.gov/sra/sra-instant/reads/ByRun/litesra/SRR/SRR036/SRR036919/SRR036919.sra 


PhiX control lane, described at:
http://www.ncbi.nlm.nih.gov/sra/SRX017204
SRR036919
We can download from here 
wget ftp://ftp.ncbi.nih.gov/sra/sra-instant/reads/ByRun/litesra/SRR/SRR036/SRR036919/SRR036919.sra 


