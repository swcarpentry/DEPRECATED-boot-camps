import animals


def test_read_animals():
    date, time, species, count = animals.read_animals('animals.txt')

    true_date = ['2011-04-22', '2011-04-23', '2011-04-23', '2011-04-23', '2011-04-23']
    true_time = ['21:06', '14:12', '10:24', '20:08', '18:46']
    true_species = ['Grizzly', 'Elk', 'Elk', 'Wolverine', 'Muskox']
    true_count = [36, 25, 26, 31, 20]

    assert date == true_date, 'Dates not the same!'
    assert time == true_time, 'Times not the same!'
    assert species == true_species, 'Species not the same!'
    assert count == true_count, 'Counts not the same!'


def test_mean():
    mean = animals.mean([4, 5, 6])
    assert mean == 5

    mean = animals.mean([10])
    assert mean == 10

    mean = animals.mean([10, 7])
    assert mean == 8.5


def test_filter_animals():
    date, time, species, count = animals.read_animals('animals.txt')
    kind = 'Grizzly'
    d, t, s, c = animals.filter_animals(kind, date, time, species, count)

    assert d == ['2011-04-22']
    assert t == ['21:06']
    assert s == ['Grizzly']
    assert c == [36]

    kind = 'Elk'
    d, t, s, c = animals.filter_animals(kind, date, time, species, count)

    assert d == ['2011-04-23', '2011-04-23']
    assert t == ['14:12', '10:24']
    assert s == ['Elk', 'Elk']
    assert c == [25, 26]


def test_mean_animals():
    mean = animals.mean_animals('animals.txt', 'Grizzly')
    assert mean == 36

    mean = animals.mean_animals('animals.txt', 'Elk')
    assert mean == 25.5
