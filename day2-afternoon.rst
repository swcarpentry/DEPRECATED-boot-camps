Shell Basics:
=====================================

A few commands used in the shell in the afternoon. 

Useful commands:
--------------------------------

- Which user is logged on. In the virtual box you should see that you are software carpentry. 

	whoami

- Unix for "print working directory". 

	pwd  

- List directory

	ls 

	ls -a  #"list all"
		
	ls -la  # long all

-Change directories:
		
	cd mydir #change directories

	cd .. # goes to directory one up from the working directory.
		
- To make a new directory:
	mkdir 2013-02-26
		
- If we want to make a new directory not nested in the one we are in we can also do that 
	mkdir /users/newstuff/

- Create a file in text editor nano:

	nano Seattle.txt

		
- Save by hitting Ctrl O and then Ctrl -X to exit. 
		
- Quick way to read file in shell
		
  cat data.txt
		
- Copying- Let's say that we need to make a copy of this file. We will type the name of the file and then its new name's. Let's make a copy and call it something more descriptive. So in this case the command is "cp." So first give it the name of the file you are trying to copy and then the new file name for the copy. 
		
	cp data.txt  mammal_abundance_UW.txt
		
- Moving files- if we want to move a file from one spot to another we can use the command mv. This is also how we rename files. 
		
	mv data.txt  /Workshop
		
	mv data.txt ..

		
- To use mv to rename we can give the destination a new filename
		
	mv data.txt data_current.txt
		
- To delete something we can use rm-for remove. And just a note that there is no "trash" or "recycling" bin in the unix shell.  Type rm file to remove a file or rmdir directory to remove a directory if it is empty. 
		
		rm data.txt
		
		rmdir emptyfolder 


Dealing with data in a file
-------------------------------------

- To sort this by column  from highest to lowest. 
	
	sort -k 2 -t, -n mammal_abundance_UW.txt
	
- Lowest to highest:
  sort -k 2 n -r mammal_abundance_UW.txt
	
- We can find out many things about these unix commands on the internet, but also by going to the manual that comes with the programs. 
	
	man sort 
	
- Exit by pressing "q".
	
	
-Save this sorted file to a new file:
	
	sort -k 1 mammal_abundance_UW.txt > sorted_mammal_abundance_UW.txt
	
- Read first 10 lines of file:
	
	head sorted_mammal_abundance_UW.txt
	
- Top 2 lines:

	head -2 sorted_mammal_abundance_UW.txt
	
- Bottom 10 lines:
	
	tail sorted_mammal_abundance_UW.txt

- Bottom 2 lines:

	tail -2 sorted_mammal_abundance_UW.txt

- Pipes: send output from one command as input to another:

	sort -k 2 n mammal_abundace_UW.txt | head -1 
	
- If we want to send the final output to a file we can also do that:
	
	sort -k 2 n mammal_abundace_UW.txt | head -1 > most_abundant_mammal_UW.txt
	
- Cut column from a file:

  cut -d, -f 1 data.txt

- Use word count wc to count the number of lines in a file
	
	wc -l  mammal_abundance_UW.txt 
	
- For multiple files:
	
	wc -l  mammal_abundance_UW.txt, mammal_abundance_UBC.txt
		
	wc -l mammal_abundance_*.txt
	
- See history of commands run in the shell session:

  history
	
- If we want we can send the print out of this command to a file. Good for our record keeping:
	
	history > command_log_2013_02_26.txt
	
- Find all lines in several files that contained "squirrels" and then aggregate these:
		grep squirrels mammal_abundance_*.txt



For loops: 
-------------------------------

	for datafile in data_*
	    do
	        echo $datafile
	        sort -k 2 -n $datafile
	    done
	 
- To make as a bash script save commands in a text editor as a .sh file and run in bash. 
	
	bash sorting_abundances.sh
	
- What about another loop with numbers 1-4?
	
	for number in {1,2,3,4}
	    do
	             echo $number
	    done 
	
- Even better would be:
	
	for number in {1..4}
	    do
	             echo $number
	    done 
	


