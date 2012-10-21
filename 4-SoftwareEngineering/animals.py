def read_animals(filename):
    """
    This function will read animal data from a file.
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
