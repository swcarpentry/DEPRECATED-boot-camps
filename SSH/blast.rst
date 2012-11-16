Installing and running standalone BLAST 

BLAST - Basic Local Alignment Search Tool


****  Download the BLAST program ****
Use an FTP client or Firefox and go to 
ftp://ftp.ncbi.nih.gov/blast/executables/

Click on 'LATEST'

For Macs you can use the .dmg file
Scroll to the bottom of the page and click on or transfer the *.dmg file
Right now that file is
ncbi-blast-2.2.25+.dmg

Go find that file and double click it
A folder will open with the ncbi-blast-2.2.25+.pkg file in it
  double click on that file and click through the steps to install it

It's installed!

Modify the .ncbirc file that was just created in your directory to have it point to the place where you'll be creating your databases

Open a terminal and at the command line type

pico .ncbirc

Make the changes
e.g.

[BLAST]
BLASTDB=/Applications/ncbi-blast-2.2.24+/db


To save and exit type
Ctrl-X
when it asks you if you want to save, type Y and then Enter when it asks if you want to save to .ncbirc



**** Creating a blast database ****


Now you can create your own BLAST databases
You can create a BLAST database from any FASTA file

There's some additional information in the README_BLAST file

There are preformatted ones for things like nr on NCBI
ftp://ftp.ncbi.nih.gov/blast/db/

and you can download any sequenced microbial genome
ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/

If the databases are really big, like they are for nr, it's not something you'll want to do on your local computer.  A server or the HPCC is better for that.


----  To create your own database from your own FASTA file ----

Copy the file into the blast database directory you just referenced in the .ncbirc

Here we'll use a test FASTA file (provided)

test_fasta.fasta

The 'formatdb' command creates the databases

-i  is the input file
-n is the name that you want for your database
-p is the type of file (protein T, or nucleotide F)

If you type 

formatdb - 
you get all the options

Run the formatdb command with both -p T and -p F so you get both the nucleotide and protein database.  The different blast programs require the different databases.

formatdb -i test_fasta.fasta -n test_fasta -p T

and

formatdb -i test_fasta.fasta -n test_fasta -p F

Now if you look in that directory you have new files and those are your blast databases


**** Running BLAST ****

The command to run BLAST is 'blastall'

If you type 
blastall - 
You get all the options

This is a standard blastall command

blastall -e 1e-05 -p tblastx -d test_fasta -i seqs.fa  -o seqs.blast

-e is the e-value cutoff you want to use.  Any matches higher than that will not be returned
-p is the program - tblastx, blastx, blastn or tblastn
-d is the database
-i is the input file
-o is the output file
-m is the output type you want
   If you're parsing the output, then you want to use -m 8.  It outputs a tab delimited format that's easy to look through
   The default shows you all the alignments

If you do use -m 8 this is the information in each column

# Query id # Subject id # % identity # alignment length # number of mismatches # number of gap openings # position of query start # position of query end # position of subject start # position of subject end # e-value of a hit # bit score of a hit  


That's it, now you have your blast information and you can parse the BLAST output
