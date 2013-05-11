import nose.tools as nt
from meananimals import read_file, calc_mean, filter_animals, mean_animals

def test_read_animals():
    date, time, animal, count = read_file('animals.txt')
    ref_date = ['2011-04-22', '2011-04-23', '2011-04-23', '2011-04-23', '2011-04-23']
    ref_time = ['21:06', '14:12', '10:24', '20:08', '18:46']
    ref_animal = ['Grizzly', 'Elk', 'Elk', 'Wolverine', 'Muskox']
    ref_count = [36, 25, 26, 31, 20]
    
    assert date == ref_date, 'Dates do not match!'
    assert time == ref_time, 'Times do not match!'
    assert animal == ref_animal, 'Animals do not match!'
    assert count == ref_count, 'Counts do not match!'

def test_mean1():
    m = calc_mean([1, 2, 3])
    assert m == 2
    
def test_mean2():
    m = calc_mean([1])
    assert m == 1

def test_mean3():
    m = calc_mean([3.4, 3.5, 3.6])
    assert m == 3.5

@nt.raises(ValueError)
def test_mean4():
    m = calc_mean([])

def test_filter_animals1():
    date, time, animal, count = read_file('animals.txt')
    fdate, ftime, fanimal, fcount = filter_animals('Elk', date, time, animal, count)
    
    assert fdate == ['2011-04-23', '2011-04-23']
    assert ftime == ['14:12', '10:24']
    assert fanimal == ['Elk', 'Elk']
    assert fcount == [25, 26]

def test_mean_animals1():
    m = mean_animals('animals.txt', 'Elk')
    assert m == 25.5

def test_mean_animals2():
    m = mean_animals('animals.txt', 'Grizzly')
    assert m == 36