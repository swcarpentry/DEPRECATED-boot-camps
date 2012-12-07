
Questionnaire graphing scripts and data
=======================================

This directory holds raw data from pre- and post-boot camp questionnaires along with gnuplot scripts to turn the data into PDF graphs. 

The gnuplot scripts are embedded within shell scripts to allow command-line arguments to be passed in.

The Makefile will create PDFs in pdfs/ from the .dat files in data/. To create all the PDFs, do:

::

 make all

Thanks to Phil Fowler for help with the gnuplot scripts.

http://www.gnuplot.info/
