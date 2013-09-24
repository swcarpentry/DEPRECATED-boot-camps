def read_file(ifile):
    open_file = open(ifile, 'r')
    
    time = []
    date = []
    animal = []
    count = []
    
    for iline in open_file:
        s = iline.split()
        date.append(s[0])
        time.append(s[1])
        animal.append(s[2])
        count.append(int(s[3]))
    
    open_file.close()
    
    return date, time, animal, count

def calc_mean(ilist):
    '''
    returns the mean of input list
    '''
    if len(ilist) == 0:
        raise ValueError('Input is empty.')
        
    total = 0.0
    for num in ilist:
        total = total + num
    return total/float(len(ilist))

def filter_animals(species, date, time, animal, count):
    """
    Given a particular species, filter out the data for just that species.

    Returns four lists: date, time, animal, count.
    """
    fdate = []
    ftime = []
    fanimal = []
    fcount = []
    
    for d, t, a, c in zip(date, time, animal, count):
        if a == species:
            fdate.append(d)
            ftime.append(t)
            fanimal.append(a)
            fcount.append(c)
    
    return fdate, ftime, fanimal, fcount

def mean_animals(filename, species):
    d, t, a, c = read_file(filename)
    d, t, a, c = filter_animals(species, d, t, a, c)
    return calc_mean(c)