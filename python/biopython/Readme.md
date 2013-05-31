##How to program if you must##

If a bioinformatic task is relatively common, it is likely that **someone else has written software to do it already**.

Format conversions, paired-read merging, quality trimming, recruitment to references, diploid snp calling, haploid
re-sequencing--these are all problems that can be solved by finding out what software purports to do the job, 
figuring out how to install it, getting the data into the right input formats and getting the results out of 
whatever formats the tools write to.  

Depending on the complexity of the task and the ease-of-use and scope of the existing alternatives, it could be 
easier to adapt an existing software library or it could be easier to write the code to do it yourself.

This section describes programming using the Biopython libraries when you must.
In general you will have your own data, you will need to change its format and do stuff to it, 
you will get some reference data, perform some comparison operation, and then perform
operations on the results of the comparison in view of your hypotheses.

Solving any particular data manipulation problem takes a certain amount of time, time that should include the time to 
fix the mistakes and confirm that your code is really is doing what you think it is.   It is dramatically easier to write a program that
you will use once and then throw away than to write a proram that is useful to you again and again.  Ask: do you really
want to solve a particular (easy) format-conversion problem six times, once for each new collaboration and each new dataset?  
If you invest effort in identifying what you need and what parts you will use again and again, you can forget the details of how you solved
this (boring) problem in the first place and direct your time to more interesting things.

##An anecdote: using python to get and plot data from a web interface##
One day, a colleague of mine showed me that the MG-RAST website had an interface that would deliver a bundle of data about a dataset in response to an HTTP request.  Specifically, the request 
http://api.metagenomics.anl.gov/metagenome_statistics/mgm4440613.3?verbosity=full
has tables of numbers representing the length distribution, GC-content, and high-level summaries of the taxonomic annotations of an NGS dataset.  There's lots of good data in here.  It's encoded in JSON format. http://en.wikipedia.org/wiki/JSON 
Fortunately, there is a python module to painlessly parse JSON into a python dict of dict.
The script ```metagenome_statistics-example.py``` contains example code that retrieves data from the website, gets some of the data out of the JSON structure, and plots it.  Python code that solves the sub-problems (retrieving data via HTTP, getting data out of JSON objects, and plotting) has already been written, so I spend my time invoking and debugging calls to these subroutines instead of finding out how to write a HTTP client or a JSON parser.  

##Biopython##
Biopython has a large collection of subroutines that do potentially useful things with biological data.
Biopython is described in *Biopython: freely available Python tools for computational molecular biology and bioinformatics.* Bioinformatics. 2009 Jun 1;25(11):1422-3. doi: 10.1093/bioinformatics/btp163. http://www.ncbi.nlm.nih.gov/pubmed/19304878?dopt=Abstract
and has a detailed *Biopython Cookbook and Tutorial*
http://biopython.org/DIST/docs/tutorial/Tutorial.html describing many of the things that it does.  

First we will show an example of getting data from NCBI, show examples of how Biopython goes through data and what types of data objects it uses, and we will give some exercises for doing something to sequences one at a time.
```Bio.SeqIO.parse()``` returns a generator.  Looping through this will produce ```SeqRecord``` objects.  These have ```Seq``` and ```SeqFeature``` objects inside of them.   

###Get reference data--NCBI's EUTILS###
Sequence comparison is at the heart of bioinformatics; to do useful comparisons, you need data (sequences) against which to compare your new, exciting sequences.  NCBI provides an interface to allow automated download of various (sequence) data products using HTTP GET requests.

The documentation for this interface, called EUTILS, is here:
http://www.ncbi.nlm.nih.gov/books/NBK25499/

Pretty much anything you can do at an NCBI web site or search engine can be done using EUTILS -- retrieving pubmed abstracts, searching SRA for NGS datasets, searching for and retrieving reference genomes.  Here we will describe retrieving sequence data (protein sequences, genome sequences, or genomes with annotations) using ```efetch```.  Other routines are described in Eutils.md.

Biopython has subroutines that take EFETCH's options and parameters and return python objects.
This saves us from having to write code that directly talks to NCBI's EFETCH API, freeing us to spend our time elsewhere.  
We just need to find out how to use these subroutines.

```retrievegbk.py``` contains an example of using the ```Entrez.efetch```  biopython method to retrieve genbank records by specifying their accession number.  

```python
'''This is an example of Entrez.efetch. http://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_EFetch_ '''
#!/usr/bin/env python
from Bio import Entrez

accessionno = "NC_012920"   # This is the human mitochrondrion reference sequence 
Entrez.email = "swc@example.com        #  Always tell NCBI who you are.  And fix the SytnaxError.
Entrez.tool = "SoftwareCarpentryBootcamp"
assert Entrez.email != None     # Check that you told NCBI who you are before continuing

# Entrez.efetch is the subroutine; 
# db, id, rettype, and retmode are parameters of EFETCH
# and read() is the method that gives us the output as one big string.

genbankdata = Entrez.efetch(db="nucleotide", id=accessionno, rettype="gb", retmode="text").read()

print genbankdata

# Let's write this to a file
filename = "%s.gbk" % accessionno
open(filename, "w").write(genbankdata)
```

Note that all we did was get the data and dump it to a file here; we will go through the data and look at what is inside later.

Exercise:
```ladyslipperITSaccessionnumbers.txt``` contains 94 accession numbers for the ITS ribosomal marker sequences of certain lady slipper orchids.
Modify the recipe above to request FASTA format instead of genbank format and then find a way to download all 94 sequences.  Once you have the sequences, you can concatenate them and run your favorite multiple-sequence-alignment program.

#### Ok, just give me the data####
Easy Exercise:
Use ```retreivegbk.py``` to get NC_001422, the genome of Phi-x174, in genbank format.  We will use this below. 

### Iterating through data records ###

Biopython provides a variety of methods for stepping through data sources one record at a time.  There is variety in the places we
can get the data from, variety in the data types, and variety in the procedures used to access the data.  

+ Data sources can be web interfaces, filenames, or file handles.  

+ Data types can include amino acid sequences, short-read nucleic acid sequences with or without qualities, draft genomes in hundreds of contigs,  or complete genomes with gene coordinates, translations, and additional notes about how the genes were identified.    These are normalized into ```SeqRecord```.

+ The access procedures include opening a data source and loading it all into memory at once, or reading a data source one record at a time. 

The ```Bio.SeqIO.parse``` method returns a *generator*, an python object with methods to let us access the data.  We can put this generator in a ```for``` loop and access each record one at a time, or we can call ```list(generator)``` to load all the records into memory at once (if we have enough memory to do so).   If we don't have enough memory load all the data at once, we need to find ways to write our programs that don't require us to do so.
Note: ```SeqIO.parse``` takes the format as a mandatory second parameter.  fasta, fastq, genbank, and embl are among the supported formats.

```SeqIO.parse()``` returns ```SeqRecord``` objects, a general Biopython data type that can accomodate fasta, fastq, genbank, and other types of data, though with different levels of detail.

We can open ```data/tiny-fasta.fasta``` using ```SeqIO.parse()```, get the first record using the ```next()``` method, and find out what data type we got using ```type```:

```python
from Bio import SeqIO
firstrecord = SeqIO.parse("data/tiny-fasta.fasta", "fasta").next()
print type(firstrecord)
    Bio.SeqRecord.SeqRecord
```
As a troubleshooting mechanism, python's ```dir``` method will list the methods and atributes that are defined for this SeqRecord:
```
print dir(firstrecord)
```
We can print each of these, and find out whether each is a method or an attribute.
The attributes, ```firstrecord.id``` and ```firstrecord.seq``` contain the data we are looking for.
Other attributes such as ```firstrecord.annotations``` and ```firstrecord.features``` are only populated for input data types richer than fasta.

We can loop through all the records, one at a time, using ```for```:
```python
generator = SeqIO.parse("data/tiny-fasta.fasta", "fasta")
for seqrecord in generator:
    print seqrecord
```

This print statement gives us a human-readable summary of the contents of the seqrecord, but not necessarily the whole thing. 

While we might expect ```seqrecord.seq``` to return a string, it returns an object of the type ```Bio.Seq.Seq```.  Strings have 
methods like ```find()``` and ```upper()```.  ```Seq``` objects additionally have methods like ```reverse_complement()``` and ```translate()```.

```python
generator = SeqIO.parse("data/tiny-fasta.fasta", "fasta")
for seqrecord in generator:
    print seqrecord.id, len(seqrecord.seq)
    print str(seqrecord.seq)
    print str(seqrecord.seq.reverse_complement())
```

Finally, you can get reminders of the recipes for accessing data and descriptions of how to use SeqIO from the manual page on ```SeqIO```:
```python
help(SeqIO)
```

```SeqRecord``` have the attributes ```SeqRecord.id```, which is a string
+ ``
The oWhat?  I just iterated over ```SeqRecord``` objects and instead of getting a sequence I get a ```Seq``` object?  

Exercise:  
Modify the existing program ```exercise-reversecomplement.py``` to output fasta whose sequences have been reverse-complemented.
##How to program if you must##

If a bioinformatic task is relatively common, it is likely that **someone else has written software to do it already**.

Format conversions, paired-read merging, quality trimming, recruitment to references, diploid snp calling, haploid
re-sequencing--these are all problems that can be solved by finding out what software purports to do the job, 
figuring out how to install it, getting the data into the right input formats and getting the results out of 
whatever formats the tools write to.  

Depending on the complexity of the task and the ease-of-use and scope of the existing alternatives, it could be 
easier to adapt an existing software library or it could be easier to write the code to do it yourself.

This section describes programming using the Biopython libraries when you must.
In general you will have your own data, you will need to change its format and do stuff to it, 
you will get some reference data, perform some comparison operation, and then perform
operations on the results of the comparison in view of your hypotheses.

Solving any particular data manipulation problem takes a certain amount of time, time that should include the time to 
fix the mistakes and confirm that your code is really is doing what you think it is.   It is dramatically easier to write a program that
you will use once and then throw away than to write a proram that is useful to you again and again.  Ask: do you really
want to solve a particular (easy) format-conversion problem six times, once for each new collaboration and each new dataset?  
If you invest effort in identifying what you need and what parts you will use again and again, you can forget the details of how you solved
this (boring) problem in the first place and direct your time to more interesting things.

##An anecdote: using python to get and plot data from a web interface##
One day, a colleague of mine showed me that the MG-RAST website had an interface that would deliver a bundle of data about a dataset in response to an HTTP request.  Specifically, the request 
http://api.metagenomics.anl.gov/metagenome_statistics/mgm4440613.3?verbosity=full
has tables of numbers representing the length distribution, GC-content, and high-level summaries of the taxonomic annotations of an NGS dataset.  There's lots of good data in here.  It's encoded in JSON format. http://en.wikipedia.org/wiki/JSON 
Fortunately, there is a python module to painlessly parse JSON into a python dict of dict.
The script ```metagenome_statistics-example.py``` contains example code that retrieves data from the website, gets some of the data out of the JSON structure, and plots it.  Python code that solves the sub-problems (retrieving data via HTTP, getting data out of JSON objects, and plotting) has already been written, so I spend my time invoking and debugging calls to these subroutines instead of finding out how to write a HTTP client or a JSON parser.  

##Biopython##
Biopython has a large collection of subroutines that do potentially useful things with biological data.
Biopython is described in *Biopython: freely available Python tools for computational molecular biology and bioinformatics.* Bioinformatics. 2009 Jun 1;25(11):1422-3. doi: 10.1093/bioinformatics/btp163. http://www.ncbi.nlm.nih.gov/pubmed/19304878?dopt=Abstract
and has a detailed *Biopython Cookbook and Tutorial*
http://biopython.org/DIST/docs/tutorial/Tutorial.html describing many of the things that it does.  

First we will show an example of getting data from NCBI, show examples of how Biopython goes through data and what types of data objects it uses, and we will give some exercises for doing something to sequences one at a time.
```Bio.SeqIO.parse()``` returns a generator.  Looping through this will produce ```SeqRecord``` objects.  These have ```Seq``` and ```SeqFeature``` objects inside of them.   

###Get reference data--NCBI's EUTILS###
Sequence comparison is at the heart of bioinformatics; to do useful comparisons, you need data (sequences) against which to compare your new, exciting sequences.  NCBI provides an interface to allow automated download of various (sequence) data products using HTTP GET requests.

The documentation for this interface, called EUTILS, is here:
http://www.ncbi.nlm.nih.gov/books/NBK25499/

Pretty much anything you can do at an NCBI web site or search engine can be done using EUTILS -- retrieving pubmed abstracts, searching SRA for NGS datasets, searching for and retrieving reference genomes.  Here we will describe retrieving sequence data (protein sequences, genome sequences, or genomes with annotations) using ```efetch```.  Other routines are described in Eutils.md.

Biopython has subroutines that take EFETCH's options and parameters and return python objects.
This saves us from having to write code that directly talks to NCBI's EFETCH API, freeing us to spend our time elsewhere.  
We just need to find out how to use these subroutines.

```retrievegbk.py``` contains an example of using the ```Entrez.efetch```  biopython method to retrieve genbank records by specifying their accession number.  

```python
'''This is an example of Entrez.efetch. http://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_EFetch_ '''
#!/usr/bin/env python
from Bio import Entrez

accessionno = "NC_012920"   # This is the human mitochrondrion reference sequence 
Entrez.email = "swc@example.com        #  Always tell NCBI who you are.  And fix the SytnaxError.
Entrez.tool = "SoftwareCarpentryBootcamp"
assert Entrez.email != None     # Check that you told NCBI who you are before continuing

# Entrez.efetch is the subroutine; 
# db, id, rettype, and retmode are parameters of EFETCH
# and read() is the method that gives us the output as one big string.

genbankdata = Entrez.efetch(db="nucleotide", id=accessionno, rettype="gb", retmode="text").read()

print genbankdata

# Let's write this to a file
filename = "%s.gbk" % accessionno
open(filename, "w").write(genbankdata)
```

Note that all we did was get the data and dump it to a file here; we will go through the data and look at what is inside later.

Exercise:
```ladyslipperITSaccessionnumbers.txt``` contains 94 accession numbers for the ITS ribosomal marker sequences of certain lady slipper orchids.
Modify the recipe above to request FASTA format instead of genbank format and then find a way to download all 94 sequences.  Once you have the sequences, you can concatenate them and run your favorite multiple-sequence-alignment program.

#### Ok, just give me the data####
Easy Exercise:
Use ```retreivegbk.py``` to get NC_001422, the genome of Phi-x174, in genbank format.  We will use this below. 

### Iterating through data records ###

Biopython provides a variety of methods for stepping through data sources one record at a time.  There is variety in the places we
can get the data from, variety in the data types, and variety in the procedures used to access the data.  

+ Data sources can be web interfaces, filenames, or file handles.  

+ Data types can include amino acid sequences, short-read nucleic acid sequences with or without qualities, draft genomes in hundreds of contigs,  or complete genomes with gene coordinates, translations, and additional notes about how the genes were identified.    These are normalized into ```SeqRecord```.

+ The access procedures include opening a data source and loading it all into memory at once, or reading a data source one record at a time. 

The ```Bio.SeqIO.parse``` method returns a *generator*, an python object with methods to let us access the data.  We can put this generator in a ```for``` loop and access each record one at a time, or we can call ```list(generator)``` to load all the records into memory at once (if we have enough memory to do so).   If we don't have enough memory load all the data at once, we need to find ways to write our programs that don't require us to do so.
Note: ```SeqIO.parse``` takes the format as a mandatory second parameter.  fasta, fastq, genbank, and embl are among the supported formats.

```SeqIO.parse()``` returns ```SeqRecord``` objects, a general Biopython data type that can accomodate fasta, fastq, genbank, and other types of data, though with different levels of detail.

We can open ```data/tiny-fasta.fasta``` using ```SeqIO.parse()```, get the first record using the ```next()``` method, and find out what data type we got using ```type```:

```python
from Bio import SeqIO
firstrecord = SeqIO.parse("data/tiny-fasta.fasta", "fasta").next()
print type(firstrecord)
    Bio.SeqRecord.SeqRecord
```
As a troubleshooting mechanism, python's ```dir``` method will list the methods and atributes that are defined for this SeqRecord:
```
print dir(firstrecord)
```
We can print each of these, and find out whether each is a method or an attribute.
The attributes, ```firstrecord.id``` and ```firstrecord.seq``` contain the data we are looking for.
Other attributes such as ```firstrecord.annotations``` and ```firstrecord.features``` are only populated for input data types richer than fasta.

We can loop through all the records, one at a time, using ```for```:
```python
generator = SeqIO.parse("data/tiny-fasta.fasta", "fasta")
for seqrecord in generator:
    print seqrecord
```

This print statement gives us a human-readable summary of the contents of the seqrecord, but not necessarily the whole thing. 

While we might expect ```seqrecord.seq``` to return a string, it returns an object of the type ```Bio.Seq.Seq```.  Strings have 
methods like ```find()``` and ```upper()```.  ```Seq``` objects additionally have methods like ```reverse_complement()``` and ```translate()```.

```python
generator = SeqIO.parse("data/tiny-fasta.fasta", "fasta")
for seqrecord in generator:
    print seqrecord.id, len(seqrecord.seq)
    print str(seqrecord.seq)
    print str(seqrecord.seq.reverse_complement())
```

Finally, you can get reminders of the recipes for accessing data and descriptions of how to use SeqIO from the manual page on ```SeqIO```:
```python
help(SeqIO)
```

```SeqRecord``` have the attributes ```SeqRecord.id```, which is a string
+ ``
The oWhat?  I just iterated over ```SeqRecord``` objects and instead of getting a sequence I get a ```Seq``` object?  

Exercise:  
Modify the existing program ```exercise-reversecomplement.py``` to output fasta whose sequences have been reverse-complemented.

####Fasta sequence parsing####
The minimal data type for sequence data, this format includes only a text record description and a (possibly long) sequence.  Nucleic acid sequences, partially-ambiguous nucleic acid sequences, and amino acid sequences can all be encoded in this bare-bones format.  

```SeqRecord``` data types have the attributes
+ ```.name```  which is the **fasta id** -- all the text before the first whitespace on the header line
+ ```.description``` the entire header line, including the fasta id and anything after the first whitespace
+ ```.seq```  the sequence, as a ```SeqIO.Seq``` object

```python
from Bio import SeqIO
for seqrecord in SeqIO.parse("tiny-fasta.fasta", "fasta"):
     print seqrecord.name, seqrecord.seq
```

Exercise:  
Modify the existing program ```exercise-allsixframes.py``` to read in nucleic acid FASTA and output six sequences, representing the translation of each sequence in all six reading frames.  Hint:  Slicing works on ```Seq``` objects like it does on strings, so if ```seq``` is one of the ```seqrecord.seq``` objects, 
```seq[1:]``` and ```seq[2:]``` the sequence with the first character chopped off, and the sequence with the first two characters chopped off, respectively.  
Hint: `Seq` objects also have  `seq.reverse_complement()` and `seq.translate()` methods that return ```Seq``` objects with the reverse complement and translation, defaulting to the standard prokaryotic code.

If you want to manipulate (say, output) the ```Seq``` objects yourself, ```str(seqrecord.seq)``` will return a string representation of the sequence.  

####Fastq sequence parsing####
FASTQ is the de-facto standard data format for the output of modern high-throughput sequencing machines.
In addition to sequence ID and header, this format includes a quality symbol for each base.

Fastq can be parsed using exactly the same ```SeqIO.parse()``` method, just with ```fastq``` instead of ```fasta``` as the format parameter.
The sequence is in the ```seq``` attribute and the quality scores (as a list of ints) is in the ```letter_annotations["phred_quality"]``` attribute.    This snippet just loops through a small example fastq file and prints the data fields:
```python
from Bio import SeqIO
generator = SeqIO.parse("data/tiny-fastq.fastq", "fastq")
for sequence in generator:
     print sequence.id
     print sequence.description
     print sequence.seq
     print sequence.letter_annotations["phred_quality"]
```
This snippet will shorten the sequences to include only the first 30 base pairs of each read:
```python
from Bio import SeqIO
generator = SeqIO.parse("data/tiny-fastq.fastq", "fastq")
for sequence in generator:
     shortsequence = sequence[0:30]
     sys.stdout.write(shortsequence.format("fastq"))
```

Exercise:  
`exercise-b-trim.py` contains a fastq parser.  Write a subroutine to perform "B-trimming" -- removing bases that have very low quality scores from the end of the reads.  You will need to determine how many bases at the end of the read have quality scores of 2 or below and remove them.   The exercise script has the input and the output in place.

####Genbank sequence parsing####

Unlike FASTA and FASTQ, the Genbank format has many optional fields, and the optional fields may be of variable lengths.  
Consequently, the `SeqRecord` data structures created by Biopython's genbank parser have more fields defined, including variable-length data types for things like gene annotations.

Using `NC_001422.1.gbk` as an example, let us examine the data structures that we get from parsing it.

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
```gbrecord.features[2]``` is the first protein-coding annotation "CDS".  Examining it, we see that its ```qualifiers``` attribute has most of the good stuff:
```python
firstcdsfeature = gbrecord.features[2]
print firstcdsfeature
print firstcdsfeature.qualifiers
print firstcdsfeature.qualifiers.keys()
    ['function', 'locus_tag', 'codon_start', 'product', 'transl_table', 'note', 'db_xref', 'translation', 'gene', 'protein_id']
```
We can output a table of all the CDS features, then, by looping over all the features, testing for CDS, and printing the fields we like:

```python
for i in gbrecord.features:
    if i.type == "CDS" :
        print i.qualifiers["locus_tag"][0], i.qualifiers["protein_id"][0] , i.qualifiers["gene"][0],  i.qualifiers["translation"][0] 

```
Oops.  This doesn't work. This is because not all the CDS records have ```qualifiers["gene"]``` defined, so python raises a KeyError when I try to access qualifiers["gene"] for the feature that doesn't have the key "gene".  We can correct this problem by putting the line that tries to access ```qualifiers["gene"]``` inside of a ```try... except``` conditional.  This lets us fill it in with a default value if "gene" isn't defined.
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
 

###High-throughput data--getting it###
High-throughput sequencing datasets range in size from a few megabytes to a few hundreds of gigabytes in size.  
Some institutions make raw sequence data available by FTP, but the sequence archive is the largest sequence data warehouse.

The NCBI offers a guide to downloading data here http://www.ncbi.nlm.nih.gov/books/NBK47540/
which includes links to downloading the *SRA toolkit*: http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=std

The sequence read archive maintains its own formats, and its own libary of programs to get data out of the SRA format.  The options for the utilities (and the formats themselves) change from time to time, so if something doesn't work, the first thing the help desk will ask you to do is update your copy of the sra toolkit.

To illustrate getting short-read sequencing data form SRA, let's get an illumina sequencing dataset for with the PhiX control genome described at http://www.ncbi.nlm.nih.gov/sra/SRX017204 .  The SRR accession number is SRR036919, and it's a 1x45 bp sequencing run.

We can download the SRA-formatted dataset from here 
wget ftp://ftp.ncbi.nih.gov/sra/sra-instant/reads/ByRun/litesra/SRR/SRR036/SRR036919/SRR036919.sra 
This is only a 300 Mbyte download.  You can expect NGS datasets (particularly shotgun and metatranscriptomic datasets) to be larger.  

Now we ne

###Being smart### 

Life is short and we have better things to do than solve easy problems.   

Here are some meta-strategies:

+ Test the accuracy of the procedure on data with known correct answers.   It's the only way you will know.

+ Test that all the parts work with each other first with a *small* subset of the data.  If something isn't working, you want to know now, not after ten hours.  Do as much testing and debugging as you can, when it's cheap, before scaling up the the whole zottobyte dataset.
 
+ Estimate how long your tasks are going to take, if you can.    For much of seuqnece analysis, ten times as much data take ten times as long to move, process.

+ Plan like you're going to have to do <any particular task> again.  A lot.  You probably are.  You probably made a mistake.  


####Fasta sequence parsing####
The minimal data type for sequence data, this format includes only a text record description and a (possibly long) sequence.  Nucleic acid sequences, partially-ambiguous nucleic acid sequences, and amino acid sequences can all be encoded in this bare-bones format.  

```SeqRecord``` data types have the attributes
+ ```.name```  which is the **fasta id** -- all the text before the first whitespace on the header line
+ ```.description``` the entire header line, including the fasta id and anything after the first whitespace
+ ```.seq```  the sequence, as a ```SeqIO.Seq``` object

```python
from Bio import SeqIO
for seqrecord in SeqIO.parse("tiny-fasta.fasta", "fasta"):
     print seqrecord.name, seqrecord.seq
```

Exercise:  
Modify the existing program ```exercise_allsixframes.py``` to read in nucleic acid FASTA and output six sequences, representing the translation of each sequence in all six reading frames.  Hint:  Slicing works on ```Seq``` objects like it does on strings, so if ```seq``` is one of the ```seqrecord.seq``` objects, 
```seq[1:]``` and ```seq[2:]``` the sequence with the first character chopped off, and the sequence with the first two characters chopped off, respectively.  
Hint: `Seq` objects also have  `seq.reverse_complement()` and `seq.translate()` methods that return ```Seq``` objects with the reverse complement and translation, defaulting to the standard prokaryotic code.

If you want to manipulate (say, output) the ```Seq``` objects yourself, ```str(seqrecord.seq)``` will return a string representation of the sequence.  

####Fastq sequence parsing####
FASTQ is the de-facto standard data format for the output of modern high-throughput sequencing machines.
In addition to sequence ID and header, this format includes a quality symbol for each base.

Fastq can be parsed using exactly the same ```SeqIO.parse()``` method, just with ```fastq``` instead of ```fasta``` as the format parameter.
The sequence is in the ```seq``` attribute and the quality scores (as a list of ints) is in the ```letter_annotations["phred_quality"]``` attribute.    This snippet just loops through a small example fastq file and prints the data fields:
```python
from Bio import SeqIO
generator = SeqIO.parse("data/tiny-fastq.fastq", "fastq")
for sequence in generator:
     print sequence.id
     print sequence.description
     print sequence.seq
     print sequence.letter_annotations["phred_quality"]
```
This snippet will shorten the sequences to include only the first 30 base pairs of each read:
```python
from Bio import SeqIO
generator = SeqIO.parse("data/tiny-fastq.fastq", "fastq")
for sequence in generator:
     shortsequence = sequence[0:30]
     sys.stdout.write(shortsequence.format("fastq"))
```

Exercise:  
`exercise-b-trim.py` contains a fastq parser.  Write a subroutine to perform "B-trimming" -- removing bases that have very low quality scores from the end of the reads.  You will need to determine how many bases at the end of the read have quality scores of 2 or below and remove them.   The exercise script has the input and the output in place.

####Genbank sequence parsing####

Unlike FASTA and FASTQ, the Genbank format has many optional fields, and the optional fields may be of variable lengths.  
Consequently, the `SeqRecord` data structures created by Biopython's genbank parser have more fields defined, including variable-length data types for things like gene annotations.

Using `NC_001422.1.gbk` as an example, let us examine the data structures that we get from parsing it.

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
```gbrecord.features[2]``` is the first protein-coding annotation "CDS".  Examining it, we see that its ```qualifiers``` attribute has most of the good stuff:
```python
firstcdsfeature = gbrecord.features[2]
print firstcdsfeature
print firstcdsfeature.qualifiers
print firstcdsfeature.qualifiers.keys()
    ['function', 'locus_tag', 'codon_start', 'product', 'transl_table', 'note', 'db_xref', 'translation', 'gene', 'protein_id']
```
We can output a table of all the CDS features, then, by looping over all the features, testing for CDS, and printing the fields we like:

```python
for i in gbrecord.features:
    if i.type == "CDS" :
        print i.qualifiers["locus_tag"][0], i.qualifiers["protein_id"][0] , i.qualifiers["gene"][0],  i.qualifiers["translation"][0] 

```
Oops.  This doesn't work. This is because not all the CDS records have ```qualifiers["gene"]``` defined, so python raises a KeyError when I try to access qualifiers["gene"] for the feature that doesn't have the key "gene".  We can correct this problem by putting the line that tries to access ```qualifiers["gene"]``` inside of a ```try... except``` conditional.  This lets us fill it in with a default value if "gene" isn't defined.
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
 

###High-throughput data--getting it###
High-throughput sequencing datasets range in size from a few megabytes to a few hundreds of gigabytes in size.  
Some institutions make raw sequence data available by FTP, but the sequence archive is the largest sequence data warehouse.

The NCBI offers a guide to downloading data here http://www.ncbi.nlm.nih.gov/books/NBK47540/
which includes links to downloading the *SRA toolkit*: http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=std

The sequence read archive maintains its own formats, and its own libary of programs to get data out of the SRA format.  The options for the utilities (and the formats themselves) change from time to time, so if something doesn't work, the first thing the help desk will ask you to do is update your copy of the sra toolkit.

To illustrate getting short-read sequencing data form SRA, let's get an illumina sequencing dataset for with the PhiX control genome described at http://www.ncbi.nlm.nih.gov/sra/SRX017204 .  The SRR accession number is SRR036919, and it's a 1x45 bp sequencing run.

We can download the SRA-formatted dataset from here 
wget ftp://ftp.ncbi.nih.gov/sra/sra-instant/reads/ByRun/litesra/SRR/SRR036/SRR036919/SRR036919.sra 
This is only a 300 Mbyte download.  You can expect NGS datasets (particularly shotgun and metatranscriptomic datasets) to be larger.  

Now we ne

###Being smart### 

Life is short and we have better things to do than solve easy problems.   

Here are some meta-strategies:

+ Test the accuracy of the procedure on data with known correct answers.   It's the only way you will know.

+ Test that all the parts work with each other first with a *small* subset of the data.  If something isn't working, you want to know now, not after ten hours.  Do as much testing and debugging as you can, when it's cheap, before scaling up the the whole zottobyte dataset.
 
+ Estimate how long your tasks are going to take, if you can.    For much of seuqnece analysis, ten times as much data take ten times as long to move, process.

+ Plan like you're going to have to do <any particular task> again.  A lot.  You probably are.  You probably made a mistake.  

