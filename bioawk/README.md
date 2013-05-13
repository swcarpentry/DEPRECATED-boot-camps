# A Quick bioawk tutorial

There was some interest in bioawk, a useful awk fork for handling
bioinformatics formats at the UC Davis Software Carpentry course, so
here is a quick tutorial.

## Concepts

Don't write your own FASTA/FASTQ parsers! FASTA is much easier, but
*code reuse* is important here. FASTQ is a very hard format to parse
*safely* and *quickly*. See
[this post](http://www.biostars.org/p/10353/#11256) for more info on
how tricky it can be to parse FASTQ.


Heng Li (author of samtools, bwa) has written a nice set of parsers
for different languages called
[readfq](https://github.com/lh3/readfq). bioawk is also a tool of Heng
Li's too. 

## Installing 

For this quick tutorial, let's just clone bioawk and install it to
your `/usr/local/bin/`:

    git clone git://github.com/lh3/bioawk.git && cd bioawk && make && mv awk bioawk && sudo cp bioawk /usr/local/bin/


## What's awk?

awk is an old programming language that is useful for processing
*streams* of data. Very quickly, there are three parts of each quick
awk program: BEGIN, middle, and END. BEGIN is loaded before the data,
the middle is run on each line, and the END is run after the lines are
done processing. For example, we could do:

    $ cat animals.txt
    2011-04-22 21:06 Grizzly 36
    2011-04-23 14:12 Elk 25
    2011-04-23 10:24 Elk 26
    2011-04-23 20:08 Wolverine 31
    2011-04-23 18:46 Muskox 20

    $ awk '{ if ($4 > 26) print $3 }' animals.txt
    Grizzly
    Wolverine

So awk maps these columns to column numbers, like `$1` and `$3`. `$0`
is the whole line. Delimiters (field separators) can be set with `-F`.

## Bioawk

Bioawk is just like awk, but instead of working with mapping columns
to variables for you, it maps bioinformatics field formats (like
FASTA/FASTQ name and sequence).

You can count sequences very effectively with bioawk, because awk
updates the built-in variable `NR` (number of records):

    bioawk -cfastx 'END{print NR}' test.fastq

But this is just the beginning; what if you wanted to use it to make a
tab-delimited table of names and sequence lengths, you could do:

    bioawk -cfastx '{print $name, length($seq)}' test-trimmed.fastq

Or maybe you want to see how many sequences are shorter (less than
80bp) now?

    bioawk -cfastx 'BEGIN{ shorter = 0} {if (length($seq) < 80) shorter += 1} END {print "shorter sequences", shorter}' test-trimmed.fastq
	
bioawk can also take other input formats: 

    bed:
         1:chrom 2:start 3:end 4:name 5:score 6:strand 7:thickstart 8:thickend 9:rgb 10:blockcount 11:blocksizes 12:blockstarts
    sam:
        1:qname 2:flag 3:rname 4:pos 5:mapq 6:cigar 7:rnext 8:pnext 9:tlen 10:seq 11:qual
    vcf:
        1:chrom 2:pos 3:id 4:ref 5:alt 6:qual 7:filter 8:info
    gff:
        1:seqname 2:source 3:feature 4:start 5:end 6:score 7:filter 8:strand 9:group 10:attribute
    fastx:
    	1:name 2:seq 3:qual 4:comment

## Other Resources

This is just a quick tutorial; see other resources like:

 - [bioawk Github repository](https://github.com/lh3/bioawk)
 
 - [another bioawk tutorial]http://kevin-gattaca.blogspot.com/2012/07/fwd-bioawk-awk-for-gziped-bed-gff-sam.html)
 
 
