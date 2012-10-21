def read_animals(filename):
    """
    Return four lists of the dates, times, animals, and counts
    from filename.
    
    """
    date = []
    time = []
    animal = []
    number = []
    f = open(filename, 'r')
    for line in f:
        d, t, a, n = line.split()
        date.append(d)
        time.append(t)
        animal.append(a)
        number.append(int(n))
    f.close()
    return date, time, animal, number
    

def mean(l):
    """
    Return the mean of a list of numbers.
    
    """
    if len(l) == 0:
        return 0
    
    return float(sum(l)) / len(l)

def filter_animals(date, time, animal, number, species):
    """
    Return the subset of date, time, andimal, and number that apply
    to species.

    """
    sub_date = []
    sub_time = []
    sub_animal = []
    sub_number = []

    for i, a in enumerate(animal):
        if a == species:
            sub_date.append(date[i])
            sub_time.append(time[i])
            sub_animal.append(a)
            sub_number.append(number[i])
    
    return sub_date, sub_time, sub_animal, sub_number

    
def mean_animals_sighted(filename, species):
    date, time, animal, number = \
        read_animals(filename)
    date, time, animal, number = \
        filter_animals(date, time, animal, number, species)
    
    if len(number) == 0:
        raise ValueError(species + ' is not in file.')
    
    mean_sighted = mean(number)
    return mean_sighted
    
    