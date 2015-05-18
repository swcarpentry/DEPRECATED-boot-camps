import sys
import os
import us
import math
from random import choice, randint, random
import calendar
import numpy as np

class Person:
    ## return int array of normally-distributed values
    def intNormal(mean, sd, size = 300):
        return [int(abs(round(x))) for x in np.random.normal(mean, sd, size)]

    ## random normal based on ACS data (http://ssc.wisc.edu/~ahanna/ACS2005-Codebook.pdf)
    eds     = intNormal(8.64, 3.95)
    incomes = intNormal(31807.47, 42190.11)
    ages    = intNormal(38.55, 23.01)
    hours   = intNormal(38.88, 13.06)

    states  = [x.name for x in us.states.STATES]

    ## remove ages lower than 18 and over 90
    for a in ages[:]:
        if a < 18 or a > 90:
            ages.remove(a)

    ## round income up to the closest hundred
    incomes = [int(round(x/100.)*100) for x in incomes]

    ## educational attainment must be at least 1
    for a in eds[:]:
        if a < 1:
            eds.remove(a)

    ## ensure unique ids
    serialNum = 101

    ## about half and half with one other (declined, identifies as neither)
    genders = ['M']*49 + ['F']*50 + ['O']

    def __init__(self):
        self.subject = Person.serialNum
        self.age     = choice(Person.ages)
        self.ed      = choice(Person.eds)
        self.income  = choice(Person.incomes)
        self.gender  = choice(Person.genders)
        self.hours   = choice(Person.hours)
        self.state   = choice(Person.states) 

        ## give unique ID
        Person.serialNum += 1

class Measurement:
    #incompleteFraction = 0.05
    def randomDate(self):
        hrs = range(8,17)
        mins = range(1,60)
        secs = range(1,60)
        months = range(5,10)

        month = choice(months)
        monthname = calendar.month_abbr[month]
        day = choice(range(1,calendar.monthrange(2015, month)[1]))
        dayname = calendar.day_abbr[calendar.weekday(2015, month, day)]
        hr = choice(hrs)
        min = choice(mins)
        sec = choice(secs)
        
        datestring = '%s %s %d %02d:%02d:%02d %s' % (dayname, monthname, day, hr, min, sec, '2015')
        return [datestring, month, day, hr, min, sec]

    def limit(self,n):
        if n < 1 :
            n = 1
        if n > 10 :
            n = 10
        return n 

    def __init__(self, p):
        """Generate a result"""
        self.person    = p
        self.serialNum = p.subject
        self.datestring, self.month, self.day, self.hr, self.min, self.sec = self.randomDate();

    def __str__(self):
        text = '# ' + '\n'
        text += "%s: %s\n" % ( 'Reported', self.datestring )
        text += "%s: %s\n" % ( 'ID',  self.person.subject )
        text += "%s: %d\n" % ( 'Age', self.person.age )
        text += "%s: %s\n" % ( 'Gender', self.person.gender )
        text += "%s: %s\n" % ( 'State of residence', self.person.state )
        text += "%s: %d\n" % ( 'Income', self.person.income )
        text += "%s: %d\n" % ( 'Education', self.person.ed )
        text += "%s: %d\n" % ( 'Hours per week', self.person.hours )
    
        return text

class Datataker:
    names = ['angela', 'JamesD', 'jamesm', 'Frank_Richard',\
        'lab183','THOMAS','alexander','Beth','Lawrence',\
        'Toni', 'gerdal', 'Bert', 'Ernie', 'olivia', 'Leandra',\
        'sonya_p', 'h_jackson'] 
    filenamestyles = ['data_%d','Data%04d','%d','%04d','surveyresult-%05d']
    suffixstyles = ['.dat','.txt','','','.DATA']

    def __init__(self):
        self.name = choice(Datataker.names)
        Datataker.names.remove(self.name)
        self.filenameprefix = choice(Datataker.filenamestyles)
        self.filenamesuffix = choice(Datataker.suffixstyles)
        self.measures = []

    def addmeasurement(self,measurement):
        self.measures.append(measurement)

    def write(self):
        os.mkdir(self.name)
        os.chdir(self.name)

        for m in self.measures:
            fname = self.filenameprefix % m.serialNum + self.filenamesuffix
            file = open(fname, 'w')
            file.write(str(m))
            file.close()
        os.chdir('..')
            
 
def main():
    npeople = 300 # should generate ~ .9*300 + 3.5*.1*300 ~ 375 files
    nfiles = 351

    people = []
    for pnum in range(npeople):
        people.append(Person())

    measurements = []
    for p in people:
        measurements.append(Measurement(p))

    nexperimenters = 7
    experimenters = []
    for i in range(nexperimenters):
        experimenters.append(Datataker())

    for fnum in xrange(min(len(measurements), nfiles)):
        ex = choice(experimenters)
        ex.addmeasurement(measurements[fnum]) 

    os.mkdir('data')
    os.chdir('data')
    for ex in experimenters:
        ex.write()
    os.chdir('..')

if __name__=='__main__':
    sys.exit(main())

