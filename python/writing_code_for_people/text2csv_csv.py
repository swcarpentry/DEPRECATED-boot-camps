#!/usr/bin/python 

def extractData(line):
    """Assume the first colon (:) separates the key from the value."""
    separator = ':'
    line_data = line.strip().split(separator)
    key = line_data[0]
    # the value may contain a separator so will need to be reassembled from
    # possibly multiple elements
    value = separator.join(line_data[1:]).strip()
    numeric_columns = ("CI type","Volume","Range","Discrimination")
    if key in numeric_columns and value != "" :
        value = float(value)
    return key,value


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

import sys
import csv
filelist = sys.argv[1:]
column_labels = ("Subject","Reported","Year/month of birth",
                 "Sex","CI type","Volume","Range","Discrimination")
all_data = []

for filename in filelist:
    data_record = parseFile(filename)
    all_data.append(data_record)

csv_writer = csv.DictWriter(sys.stdout, delimiter=',', fieldnames=column_labels, quoting=csv.QUOTE_MINIMAL)

csv_writer.writeheader()
csv_writer.writerows(all_data)

# for data_record in all_data:
#     csv_writer.writerow(data_record)



        


