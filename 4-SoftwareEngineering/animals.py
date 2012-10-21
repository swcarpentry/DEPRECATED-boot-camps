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


def filter_animals(date, time, animals, count, species):
    """
    Return only the data that apply to species.

    """
    filtered_date = []
    filtered_time = []
    filtered_animals = []
    filtered_count = []

    for i, a in enumerate(animals):
        if a == species:
            filtered_date.append(date[i])
            filtered_time.append(time[i])
            filtered_animals.append(animals[i])
            filtered_count.append(count[i])

    return filtered_date, filtered_time, filtered_animals, filtered_count


def mean_animals_sighted(filename, species):
    """
    Return the mean number of animals seen per sighting
    in filename.

    """
    date, time, animals, count = read_animals(filename)
    date, time, animals, count = \
        filter_animals(date, time, animals, count, species)

    if len(count) == 0:
        raise ValueError('That animal is not in the file!')

    mean_sighted = mean(count)
    return mean_sighted
