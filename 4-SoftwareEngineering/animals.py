def read_animals(filename):
    """
    This function will read animal data from a file. It returns 4 lists,
    one for each column in the file.

    """

    date = []
    time = []
    animal = []
    number = []

    f = open(filename, 'r')

    # Go through each line and append each column to the appropriate list
    for line in f:
        d, t, a, n = line.split()
        date.append(d)
        time.append(t)
        animal.append(a)
        number.append(int(n))

    f.close()
    return (date, time, animal, number)


def mean(l):
    """
    Return the mean of a list of numbers. Returned value
    should always be a float.

    """
    if len(l) == 0:
        return None

    return float(sum(l)) / len(l)
