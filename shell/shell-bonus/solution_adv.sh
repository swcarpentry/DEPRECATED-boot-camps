#!/bin/bash

#In this example we'll make a temporary file (I normally call them tmp or temp, sometimes if I'm really lazy t)

for f in `ls athletes/*.txt`
do
    # I've used the long names for the options because I feel they're much more readable
    # and I can often times better remember the name of the option rather than the
    # cursory related single character

    # We want to sort by the total number of medals won, not golds for this case
    # The $'\t' is because if bash sees the \t it'll try to expand it, this forces
    # bash to pass the tab character to sort
    sort --reverse --numeric --key=10 --field-separator=$'\t' $f | head -n 1 >> tmp_solution
done

# I want to write the output to the screen for the user
cat tmp_solution
cut -f10 tmp_solution | ./mean.py # Cut the gold medals column out of the combined results and pipe the numbers in to the mean.py program we wrote
#rm tmp_solution # clean up after ourselves