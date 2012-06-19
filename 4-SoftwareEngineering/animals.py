def read_animals(filename):
    """
    Reads animal count file agreed upon at May 30, meeting.
    Must have 4 columns, and no header row.
    Columns are date, time, animal name, number seen.
    """
    f = open(filename, 'r')
    date = []
    time = []
    animal = []
    number = []
    
    # iterate over the file one line at a time
    for line in f:
        d, t, a, n = line.split()
        date.append(d)
        time.append(t)
        animal.append(a)
        number.append(int(n))
    return date, time, animal, number

def mean(l):
    """
    Returns the mean value of the given list
    """
    sum = 0
    for x in l:
        sum = sum + x
    return sum / float(len(l))

def get_animal(date, time, animal, number, animal_name):
    """
    Given lists of dates, times, animals and numbers, return
    only those pertaining to the given animal_name.
    """
    new_date = []
    new_time = []
    new_number = []
    
    # for i in range(len(animal)):
    
    for d, t, a, n in zip(date, time, animal, number):
        if a == animal_name:
	    new_date.append(d)
	    new_time.append(t)
	    new_number.append(n)
	    
    return new_date, new_time, new_number
    
def main(filename, animal):
    dates, times, animals, counts = read_animals(filename)
    dates, times, counts = get_animal(dates, times, animals, counts, animal)
    mean_count = mean(counts)
    return mean_count
    
    