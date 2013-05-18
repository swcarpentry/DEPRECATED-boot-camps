The boot camp at NESCent May 16 - 17 incorporated lots of R content. Here's a guide to the files, how they evolved, what purposes they served.

block\_01: Students fired up RStudio. Basic exploration of the environment RStudio provides. Entered commands live in the R console. Discussed organizing a project, especially an R analytical project, and R's notions of workspace and working directory. Created an RStudio project to use for the remainder of the bootcamp. Sent commands from the History to the source editor. Eventually saved these as the stand-alone script [block01\_toyExample.R](https://github.com/jennybc/boot-camps/blob/2013-05-nescent/R/code/block01_toyExample.R). Practiced grooming code and using RStudio's facilities for sending code from the source editor to the console.
  * inputs: none
  * code: [block01\_toyExample.R](https://github.com/jennybc/boot-camps/blob/2013-05-nescent/R/code/block01_toyExample.R)
  * outputs: [avgX.txt](https://github.com/jennybc/boot-camps/blob/2013-05-nescent/R/results/avgX.txt), [niftyPlot.pdf](https://github.com/jennybc/boot-camps/blob/2013-05-nescent/R/figs/niftyPlot.pdf)
  
block\_02: Basic care and feeding of the most common R objects. Special emphasis on `data.frames`. `read.table` and friends for import. Bit of figure-making with the `lattice` package. Using the `subset()` and `with()` functions and the `data=` and `subset=` arguments found in many functions to do computations _in situ_ with added bonus of readable code. How to access various bits of various R objects, i.e. indexing. Accurate transcript can be found in the script [block02\_careFeedingData.R](https://github.com/jennybc/boot-camps/blob/2013-05-nescent/R/code/block02_careFeedingData.R), which is NOT meant to be run as a whole -- it's for interactive use.
  * inputs: [gapminderDataFiveYear.txt](https://github.com/jennybc/boot-camps/blob/2013-05-nescent/R/data/gapminderDataFiveYear.txt)
  * code: [block02\_careFeedingData.R](https://github.com/jennybc/boot-camps/blob/2013-05-nescent/R/code/block02_careFeedingData.R)
  * output: none
  
block03: Data aggregation = doing something repetitive for various logical bits of an R object. E.g. taking means of rows in a matrix, computing statistical summaries for variables in a `data.frame`, fitting a model to sub-`data.frames` induced by separating the Gapminder data out by country. Used the `apply` family of functions in base R and also introduced the add-on package `plyr`. Accurate transcript can be found in the script [block03\_dataAggregation.R](https://github.com/jennybc/boot-camps/blob/2013-05-nescent/R/code/block03_dataAggregation.R), which is NOT meant to be run as a whole -- it's for interactive use.
  * inputs: [gapminderDataFiveYear.txt](https://github.com/jennybc/boot-camps/blob/2013-05-nescent/R/data/gapminderDataFiveYear.txt)
  * code: [block03\_dataAggregation.R](https://github.com/jennybc/boot-camps/blob/2013-05-nescent/R/code/block03_dataAggregation.R)
  * output: none

block04: Sort of a capstone "putting it all together" piece. Revisiting country specific linear models of life expectancy against year. Before writing those results to file for later use, reordering the continent factor rationally (based on rate of life expectancy gains) and dropping Oceania (too few countries). Different ways to write rectangular data to file with various pros/cons: `write.table`, `dput`, `saveRDS`, which have natural relationships with `read.table`, `dget`, `readRDS`. Accurate transcript of live work can be found in the script [block04\_puttingAllTogether.R](https://github.com/jennybc/boot-camps/blob/2013-05-nescent/R/code/block04_puttingAllTogether.R), which is NOT meant to be run as a whole -- it's for interactive use. We did package some of our work nicely as scripts that could be `knit` and/or `source`'d or put into a pipeline (see below).
  * [block04\_puttingAllTogether.R](https://github.com/jennybc/boot-camps/blob/2013-05-nescent/R/code/block04_puttingAllTogether.R)
    - inputs: [gapminderDataFiveYear.txt](https://github.com/jennybc/boot-camps/blob/2013-05-nescent/R/data/gapminderDataFiveYear.txt)
    - output: none
  * [01\_countrySpecificInterceptSlope.R](https://github.com/jennybc/boot-camps/blob/2013-05-nescent/R/code/01_countrySpecificInterceptSlope.R)
    - inputs: [gapminderDataFiveYear.txt](https://github.com/jennybc/boot-camps/blob/2013-05-nescent/R/data/gapminderDataFiveYear.txt)
    - output: [gCoef.txt](https://github.com/jennybc/boot-camps/blob/2013-05-nescent/R/results/gCoef.txt), [gCoef.rds](https://github.com/jennybc/boot-camps/blob/2013-05-nescent/R/results/gCoef.rds)
  * [02\_slopeComparisonAsiaVsAmericas.R](https://github.com/jennybc/boot-camps/blob/2013-05-nescent/R/code/01_countrySpecificInterceptSlope.R) (we used this to demonstrate the super-lightweight dynamic report generation capability of RStudio: "File --> Compile notebook", also available as a button)
    - inputs: [gCoef.rds](https://github.com/jennybc/boot-camps/blob/2013-05-nescent/R/results/gCoef.rds)
    - outputs (after compiling notebook): [02_slopeComparisonAsiaVsAmericas.html](https://github.com/jennybc/boot-camps/blob/2013-05-nescent/R/prose/02_slopeComparisonAsiaVsAmericas.html)
  * [03\_slopeComparisonAsiaVsAmericas.R](https://github.com/jennybc/boot-camps/blob/2013-05-nescent/R/code/03_slopeComparisonAsiaVsAmericas.R) (we used this to demonstrate how a stand-alone script could leave files behind for later use, such as a PDF and the results of a two-sample t-test; essentially equivalent to 02\_slopeComparionsAsiaVsAmericas.R but optimized for running in a hands-off way)
      - inputs: [gCoef.rds](https://github.com/jennybc/boot-camps/blob/2013-05-nescent/R/results/gCoef.rds)
      - outputs: [slopes_AsiaVsAmericas.pdf](https://github.com/jennybc/boot-camps/blob/2013-05-nescent/R/figs/slopes_AsiaVsAmericas.pdf), [02_slopeComparisonAsiaVsAmericas_fromSink.txt](https://github.com/jennybc/boot-camps/blob/2013-05-nescent/R/prose/02_slopeComparisonAsiaVsAmericas_fromSink.txt)
      
block05: JENNY mention the stuff you'd hoped to cover on Day 1 that did not quite fit.

Jenny showed a couple slides from her UBC courses with helpful visuals for various R concepts. Students requested those and they are in `prose/slides.pdf`.

On day 2, after Ben covered command line / shell stuff, students were tasked with making subdirectories in their RStudio project directory and filing things away neatly, e.g., putting all *.R files in `code/`, and *pdf files in `/figs`. This also required updating paths for all read/write commands to reflect the new organization hierarchy.

After Karen covered some git and github stuff, students forked the nescent branch of the swcarpentry boot camps repo and got Jenny's definitive day1, post-tidy R materials.

Ben covered some `make` stuff and showed how to put much of our R work into a nice pipeline.

JENNY TO DO: understand how to achieve same as RStudio's Compile notebook from the command line, for insertion into pipelines.
