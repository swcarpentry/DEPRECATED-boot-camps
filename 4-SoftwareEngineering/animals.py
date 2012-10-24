def read_animals(filename):
    date = []
    time = []
    animal = []
    count = []
    with open(filename, 'r') as f:
        for line in f:
            d, t, a, c = line.split()
            date.append(d)
            time.append(t)
            animal.append(a)
            count.append(int(c))
    return date, time, animal, count


def mean(l):
    if not l:
        return 0

    return float(sum(l)) / len(l)


def filter_animals(date, time, animal, count, species):
    """
    Filter a set of survey data based on
    a species name. Returns four filtered lists.

    """
    r_date = []
    r_time = []
    r_animal = []
    r_count = []
    for i, a in enumerate(animal):
        if a == species:
            r_date.append(date[i])
            r_time.append(time[i])
            r_animal.append(animal[i])
            r_count.append(count[i])
    return r_date, r_time, r_animal, r_count


def mean_animals(filename, species):
    d, t, a, c = read_animals(filename)
    d, t, a, c = filter_animals(d, t, a, c, species)
    return mean(c)
