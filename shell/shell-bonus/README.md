# Automation

## Goals

* (review) turning bash commands in to shell scripts
* Using (bash) for loops to run a single command for multiple inputs
* Using make for advanced automation

## Exercises
* Write a bash script to find the line with the largest number of gold medals won in the test_data.txt file
* Modify this bash script to loop over every *.txt file in the exercise directory and find the country with the largest number of gold medals won in each file
* Use the `mean.py` script we wrote with the last bash script we wrote to compute the mean number of max medals awarded to any one country.
  * See the `cut` command for how to extract a columns from a text file.
* Use your script to compute the max medals awarded per file and the mean of the maximum

My solution to this exercise is in solution_basic.sh

## Bonus Exercises
* Make a directory `athletes`
* Using the `olympic_medals_raw.txt` file create one file per athlete
  * Suggestions:
  * Using `cut`, `sort`, and `uniq` commands, create a file containing the unique athletes in the `olympic_medals_raw.txt` file
  * Using either a bash script with `for`, `grep` and `cat` (or a python script) to output a file for each athlete
* Use your script that computes the maximum number of total medals won per year and average maximum to work with the athletes directory

My solution to this exercise is in solution_adv.sh