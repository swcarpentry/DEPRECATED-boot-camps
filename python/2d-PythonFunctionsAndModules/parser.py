f = open("data.dat")

lines = []

for line in f:
  line = line.strip()
  if line is not "":
    lines.append(line)

f.close()

dat = {}

for line in lines:
  keyval = line.split(": ")
  dat[keyval[0]] = keyval[1]
