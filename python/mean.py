import sys

def mean(ilist):
    sum = 0.0
    for num in ilist:
        sum = sum + num
    return sum / float(len(ilist))

if __name__ == "__main__":
    filename = sys.argv[1]
    ofile = open(filename, 'r')
    ilist = []
    row = ofile.read()
    for num in row.split(','):
        ilist.append(float(num))
    print mean(ilist)
