import animals


def test_read_animals_date():
    reference_date = ['2011-04-22', '2011-04-23', '2011-04-23', '2011-04-23',
                      '2011-04-23']

    date, time, animal, number = animals.read_animals('animals.txt')

    assert date == reference_date, 'Date is wrong'


def test_read_animals_time():
    reference_time = ['21:06', '14:12', '10:24', '20:08', '18:46']

    date, time, animal, number = animals.read_animals('animals.txt')

    assert time == reference_time, 'Time is wrong'


def test_mean1():
    l = [1]
    assert animals.mean(l) == 1


def test_mean2():
    l = [1, 2, 3]
    assert animals.mean(l) == 2


def test_mean3():
    l = [4, 2, 3, 5]
    assert animals.mean(l) == 3.5


def test_mean4():
    l = []
    assert animals.mean(l) == None


def test_filter_animals():
    date, time, animal, count = \
        animals.read_animals('animals.txt')

    species = 'Elk'
    date, time, animal, count = \
        animals.filter_animals(date, time, animal, count, species)

    assert count == [25, 26]
    assert animal == ['Elk', 'Elk']


def test_mean_animals_sighted1():
    mean_sighted = animals.mean_animals_sighted('animals.txt',
                                                'Wolverine')
    assert mean_sighted == 31


def test_mean_animals_sighted2():
    mean_sighted = animals.mean_animals_sighted('animals.txt',
                                                'Elk')
    assert mean_sighted == 25.5


import nose
@nose.tools.raises(ValueError)
def test_mean_animals_sighted3():
    mean_sighted = animals.mean_animals_sighted('animals.txt',
                                                'Pangolin')
