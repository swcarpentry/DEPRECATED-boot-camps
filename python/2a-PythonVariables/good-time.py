import time

print "It's (mostly) always a good time ..."

time.sleep(1)

print "to eat ..."

time.sleep(2)

f = open("cake.txt")
for line in f:
  print line.rstrip()
  time.sleep(0.2)

f.close()

