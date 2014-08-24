#!/usr/bin/python 

def extractData(line):
    """Assume the first colon (:) separates the key from the value."""
    separator = ':'
    line_data = line.strip().split(separator)
    key = line_data[0]
    # the value may contain a colon so will need to be reassembled from
    # possibly multiple elements
    value = separator.join(line_data[1:])  
    return key,value.strip()


def parseFile(filename):
    """Read all the lines from a file and return them as a dictionary of key/value pairs"""
    textfile = open(filename)
    data_record = {}
    for line in textfile:
        # all lines that contain key/value pairs must have at least on colon
        if ':' in line:
            (key,value) = extractData(line)
            data_record[key] = value
    return data_record

def writeCSVHeader(column_labels,csv_separator):
    header = []
    for column in column_labels:
        header.append('"' + column + '"')
    print csv_separator.join(header)

def writeCSVRow(column_labels,data_record,csv_separator):
    row = []
    for column in column_labels:
        row.append('"' + data_record[column] + '"')
    print csv_separator.join(row)

import sys
filelist = sys.argv[1:]
column_labels = ("Subject","Reported","Year/month of birth",
                 "Sex","CI type","Volume","Range","Discrimination")
all_data = []


for filename in filelist:
    data_record = parseFile(filename)
    all_data.append(data_record)

csv_separator = ','

writeCSVHeader(column_labels,csv_separator)

for data_record in all_data:
    writeCSVRow(column_labels,data_record,csv_separator)



        


