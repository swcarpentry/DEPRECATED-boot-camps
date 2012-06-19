import animals

def test_read_animals():
    filename = 'animals.txt'
    
    date, time, animal, count = animals.read_animals(filename)
    
    ref_date = ['2011-04-22', '2011-04-23', '2011-04-23',
                '2011-04-23', '2011-04-23']
    ref_time = ['21:06', '14:12', '10:24', '20:08', '18:46']
    ref_animals = ['Grizzly', 'Elk', 'Elk', 'Wolverine', 'Muskox']
    ref_counts = [36, 25, 26, 31, 20]
    
    assert date == ref_date, 'Dates do not match!'
    assert count == ref_counts, 'Counts do not match!'

from numpy import testing    
    
    
def test_mean():
    assert animals.mean([5]) == 5
    assert animals.mean([3, 5]) == 4
    assert animals.mean([1, 2, 3, 4]) == 2.5
    testing.assert_almost_equal(animals.mean([1, 2, 3, 4]), 2.5)

from nose.tools import raises
    
@raises(ZeroDivisionError)    
def test_mean_empty_list():   
    animals.mean([])
    
def test_get_animal():
    date, time, animal, count = animals.read_animals("animals.txt")
    date, time, count = animals.get_animal(
        date, time, animal, count, "Elk")
    
    ref_date = ['2011-04-23', '2011-04-23']
    ref_time = ['14:12', '10:24']
    ref_counts = [25, 26]
    
    assert date == ref_date, 'Dates do not match!'
    assert count == ref_counts, 'Counts do not match!'

def test_get_missing_animal():
    date, time, animal, count = animals.read_animals("animals.txt")
    date, time, count = animals.get_animal(
        date, time, animal, count, "Squirrel")
    assert date == []
    assert time == []
    assert count == []
    
def test_main():
    mean_count = animals.main("animals.txt", "Elk")
    assert mean_count == 25.5