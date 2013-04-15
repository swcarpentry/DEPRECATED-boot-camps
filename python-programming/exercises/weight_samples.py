import sys

#In this list we'll record the weights of samples.
sampleweight = []

#enterweight() prompts the user to enter weight for a sample
# and records it in the sampleweight list
def enterweight():
   weight = raw_input("Enter the sample weight (in agreed units): ")
   sampleweight.append(int(weight))
   print "Recorded sample number ",len(sampleweight)-1,"; the weigh recorded is ",  weight, " units."
   return sampleweight

#compareweight() compares weights of two selected samples
#which indexes in sampleweight[] are passed as function arguments
def compareweight(sample1,sample2):
   if sampleweight[sample1] > sampleweight[sample2]:
      print "Sample number ", sample1, " is heavier than sample number ", sample2
   elif sampleweight[sample1]==sampleweight[sample2]:
      print "Both samples are equally heavy."
   else:
      print "Sample number ", sample2," is heavier than sample number ", sample1

#sumweight() adds weights of two selected samples
#which indexes in sampleweight[] are passed as function arguments
def sumweight(sample1, sample2):
   return sampleweight[sample1] + sampleweight[sample2]

#Function enterweight() returns the list of the samples' weights (sampleweight[])
#Hence we can print the list after the function finishes.
print enterweight()
print enterweight()
print enterweight()

compareweight(0,1)
compareweight(0,2)

print "Samples ",sampleweight[0], " and ",sampleweight[1]," weigh in total ",sumweight(0,1), " units."
print 'Samples ',sampleweight[0], 'and ', sampleweight[2],' weigh in total ', sumweight(0,2), " units."
