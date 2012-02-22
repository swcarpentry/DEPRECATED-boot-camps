# find_missing.py
import os, sys, glob
import fileinput

"""These are the default values, in case an answer is missing."""
questions = {'Reported: ':'Thur Jan 1 00:00:00 1970',\
        'Subject: ':'unknown',\
        'Year/month of birth: ':'0000/00',\
        'Sex: ':'M',\
        'CI type: ':'0',\
        'Volume: ':'0',\
        'Discrimination: ':'0'\
        }

def main():
    """This method finds the incomplete tests and labels them somehow."""
    mainpath = '../cleaneddata'
    mainpattern = '*.txt'
    filelist = find_files(mainpath, mainpattern)
    for each_file in filelist : 
        file = open(each_file)
        for line in file.readlines():
            if (found_missing(line)):
                fix_missing(each_file, line)

def find_files(path, pattern):
    """This method identifies the files that contain data"""
    data = []
    filenames = glob.glob( os.path.join(path, pattern) )
    for name in filenames :
        print "current file is: " + name 
    return filenames


def found_missing(a_line):
    """This method checks that the gender line contains a missing or a Sex: N response"""
    # The questions are a dictionary of the matching string and default answer 
    to_ret = False
    for question, default in questions.items() :
        if a_line.find(question) != -1 :
            if "N" in a_line.replace(question,'').split() or len(a_line.replace(question,'').split())==0: 
                to_ret = True
                print "Found missing"
            else :
                to_ret = False
    return to_ret

def fix_missing(filename, entry):
    """This method adds the default response where one is missing"""
    for question, default in questions.items() :
        if entry.find(question) != -1 : 
            for line in fileinput.FileInput(filename, inplace=True): 
                if question == "Sex: " : 
                    line = line.replace(question + 'N', question+default)
                else : 
                    line = line.replace(question, question+default)
                sys.stdout.write(line)

class TestClass :
    """This class is a class of tests addressing the functions in this file"""

    def test_find_files(self):
        names = find_files("../data/data/jamesm", 'data_*.txt')
        assert (len(names) != 0)
        names = find_files("../data/data/jamesm", 'Data_*.txt')
        assert (len(names) == 0)
        names = find_files("fakepath", 'data_*.txt')
        assert (len(names) == 0)

    def test_not_found_missing(self):
        """This method tests the check entries function for complete lines"""
        for line in self.completelines :
            assert (not found_missing(line)) 

    def test_found_missing(self):
        """This method tests the check entries function for incomplete lines"""
        for line in self.incompletelines :
            assert (found_missing(line))

    def test_not_fix_missing(self):
        """This method tests the fix missing function for a complete file"""
        for question, default in questions.items():
            fix_missing(self.completefile.__getattribute__('name'),question)
            for line in self.completefile.readlines():
                assert (not found_missing(line))

    def test_fix_missing(self):
        """This method tests the fix missing function for an incomplete file"""
        for question, default in questions.items():
            fix_missing(self.incompletefile.__getattribute__('name'),question)
            for line in self.incompletefile.readlines():
                assert (not found_missing(line))

    def setUp(self):
        """This method sets up some fixtures"""
        # complete lines
        self.completelines = []
        for question, default in questions.items() :
            self.completelines.append(question + default)
        # incomplete lines
        self.incompletelines = []
        for question, default in questions.items() :
            self.incompletelines.append(question)
        # complete file
        self.completefile = open("testcomplete.txt",'w+')
        self.completefile.writelines(self.completelines)
        # incomplete file
        self.incompletefile = open("testincomplete.txt", 'w+ ')
        self.incompletefile.writelines(self.incompletelines)
    
    def tearDown(self):
        """This deletes the fixtures that were created in the setup"""
        os.remove(self.completefile.__getattribute__('name'))
        os.remove(self.incompletefile.__getattribute__('name'))

if __name__=='__main__':
    sys.exit(main())
