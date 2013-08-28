#!/bin/bash

for athlete in `cut -f1 olympic_athletes_raw.txt | sort | uniq`
do
    echo $athlete
    #We have to put quotes around the second athlete so that bash knows we mean the athlete variable, not athletes_medals.txt
    grep $athlete olympic_athletes_raw.txt > athletes/"$athlete"_medals.txt
done