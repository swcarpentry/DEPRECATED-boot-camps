import sys

def get_mean(ilist):
    '''
    This code calculates the mean of a list of numbers
    Input: 
        ilist: a list of integer or floating point numbers
    Output:
        mean: the mean of ilist
    '''
    total = 0
    for num in ilist:
        if (type(num) is not float) and (type(num) is not int):
            return 'Error'
        total = total + float(num)
    mean = total / float(len(ilist))
    return mean

def read_file(filename):
    '''
    This code creates a list of grades from a file with one column of student
    names and another column of student's grades
    Input: 
        filename: the name of the file containing the grades. It is assumed
        that the grades are in the second column
    Output:
        grades: a list of grades as floating point numbers
    '''
    ofile = open(filename, 'r')
    all_lines = ofile.readlines()
    grades = []
    for iline in all_lines:
        sline = iline.split()
        if len(sline) != 2:
            return 'Error'
        grades.append(float(sline[1]))
    return grades




