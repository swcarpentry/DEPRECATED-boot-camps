import animals


def test_read_animals_times():
    date, time, animal, count = \
        animals.read_animals('animals.txt')
    ref_t = ['21:06', '14:12', '10:24', '20:08', '18:46']
    assert time == ref_t, "Times don't match!"


def test_read_animals_counts():
    date, time, animal, count = \
        animals.read_animals('animals.txt')
    ref_c = [36, 25, 26, 31, 20]
    assert count == ref_c, "Counts don't match!"


def test_mean1():
    l = [1]
    assert animals.mean(l) == 1


def test_mean2():
    l = [1, 2, 3]
    assert animals.mean(l) == 2


def test_mean3():
    l = [2, 3, 5, 4]
    assert animals.mean(l) == 3.5


def test_mean4():
    l = []
    assert animals.mean(l) == 0


def test_filter_animals1():
    date, time, animal, count = \
        animals.read_animals('animals.txt')
    date, time, animal, count = \
        animals.filter_animals(date, time,
                               animal, count,
                               'Grizzly')
    assert date == ['2011-04-22']
    assert time == ['21:06']
    assert animal == ['Grizzly']
    assert count == [36]


def test_filter_animals2():
    date, time, animal, count = \
        animals.read_animals('animals.txt')
    date, time, animal, count = \
        animals.filter_animals(date, time,
                               animal, count,
                               'Elk')
    assert date == ['2011-04-23', '2011-04-23']
    assert time == ['14:12', '10:24']
    assert animal == ['Elk', 'Elk']
    assert count == [25, 26]


def test_mean_animals1():
    filename = 'animals.txt'
    species = 'Grizzly'

    assert animals.mean_animals(filename, species) == 36


def test_mean_animals2():
    filename = 'animals.txt'
    species = 'Elk'

    assert animals.mean_animals(filename, species) == 25.5
