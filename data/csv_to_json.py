"""
Copyright (c) 2013. The University of Edinburgh.

This is a simple Python script that converts a comma-separated
value file into a set of simple JSON documents that are added
to a MongoDB database.

The CSV file is assumed to be of form,

    ColumName1,ColumName2,...,ColumNameN
    Value11,Value12,...,Value1N
    Value21,Value22,...,Value2N
    ...
    ValueM1,ValueM2,...,ValueMN

Each JSON document is of form,

    {ColumnName1:ValueM1,
     ColumnName2:ValueM2,
     ...
     ColumnNameN:ValueMN}

All values are stored as strings unless they can be converted to
either floating point numbers or integers.

Usage:

    $ python csv_to_json.py CSV_FILE DATABASE_NAME COLLECTION_NAME

For example,

    $ python csv_to_json.py HYG-Database/hygxyz.csv astronomy stars

Any existing database of the same name is dropped first.
"""

import json
import sys
from pymongo import MongoClient

def get_values(line):
  """
  Given a line of comma-separated values, return the individual values
  as a list, having stripped off any leading or trailing white-space.

  @param line. Line of comma-separated values.
  @return list of values.
  """
  values = line.split(',')
  values = [value.strip() for value in values]
  return values

def get_json(columns, values):
  """
  Given a list of column names and a list of values, return a
  JSON document of form {column_1:value_1, ..., column_N:value_N}.

  @param columns: List of column names.
  @param values: List of values.
  @return JSON document.
  """
  doc = {}
  for column, value in zip(columns, values):
    try:
        value = float(value)
        if (value.is_integer()):
          value = int(value)
    except ValueError:
        pass
    doc[column] = value
  return doc

def insert_docs(file, collection):
  """
  Iterate through the rows of a comma-separated values file,
  convert each row into a JSON document and insert the document
  into a MongoDB collection. The first line of the file is
  assumed to hold column names.

  @param file: File handle.
  @param collection: MongoDB collection object.
  @return number of data lines processed.
  """
  line = file.readline()
  columns = get_values(line)
  print "Columns: ", columns
  count = 0
  # Count number of rows.
  for line in file:
    values = get_values(line)
    doc = get_json(columns, values)
    count += 1
    collection.insert(doc)
  return count

def populate_database(filename, database_name, collection_name):
  """
  Iterate through the rows of a comma-separated values file,
  convert each row into a JSON document and insert the document
  into a MongoDB collection. The first line of the file is
  assumed to hold column names.

  @param filename: File name.
  @param database_name: MongoDB database name. If the database
  already exists then it is dropped and recreated.
  @param collection_name: MongoDB collection name.
  """
  file = open(filename, 'r')
  connection = MongoClient()
  connection.drop_database(database_name)
  database = connection[database_name]
  collection = database[collection_name]
  count = insert_docs(file, collection)
  file.close()
  collection_count = collection.count()
  assert count == collection_count, "Number of documents in database does not match number of rows in file"

def print_usage():
  print "Usage: python csv_to_json.py FILE_NAME DATABASE_NAME COLLECTION_NAME"

if (len(sys.argv) < 2):
  print_usage()
  sys.exit("Missing file name")
filename = sys.argv[1]
if (len(sys.argv) < 3):
  print_usage()
  sys.exit("Missing database name")
database_name = sys.argv[2]
if (len(sys.argv) < 4):
  print_usage()
  sys.exit("Missing collection name")
collection_name = sys.argv[3]
populate_database(filename, database_name, collection_name)
