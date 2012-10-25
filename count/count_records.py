
import sys

# Given a file name, count the number of records
# in the file. Lines starting with "D" or "#"
# are ignored.
def count_records(filename):
  source = open(filename, 'r')
  count = 0
  # Count number of data records.
  for line in source:
      if line.startswith('#'): # Skip comments.
          pass        
      elif line.startswith('D'): # Skip title line.
          pass        
      else:
          count += 1
  source.close()
  return count

if (len(sys.argv) < 2):
    sys.exit("Missing file name")
filename = sys.argv[1]
print count_records(filename)
