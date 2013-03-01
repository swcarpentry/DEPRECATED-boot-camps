def mean_animals(filename, kind):
    date, time, species, count = read_animals(filename)
    date, time, species, count = filter_animals(kind, date, time, species, count)
    return mean(count)

def read_animals(filename):
    """
    Reads an animals type file. Returns four lists, one per column.

    Input is a filename.

    Order of returned items: date, time, species, count.
    """
    f = open(filename)
    dates = []
    times = []
    species = []
    counts = []

    for line in f:
        d, t, s, c = line.split()
        dates.append(d)
        times.append(t)
        species.append(s)
        counts.append(int(c))

    f.close()

    return dates, times, species, counts


def mean(nums):
    total = 0

    for n in nums:
        total = total + n

    return float(total) / len(nums)


def filter_animals(kind, date, time, species, count):
    fdate = []
    ftime = []
    fspecies = []
    fcount = []

    for i in range(len(date)):
        if species[i] == kind:
            fdate.append(date[i])
            ftime.append(time[i])
            fspecies.append(species[i])
            fcount.append(count[i])

    return fdate, ftime, fspecies, fcount
