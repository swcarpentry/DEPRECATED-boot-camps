import animals

def test_read_animals():
    date, time, animal, number = animals.read_animals('animals.txt')
  
    reference_date = ['2011-04-22', '2011-04-23', '2011-04-23', '2011-04-23', '2011-04-23']
    reference_time = ['21:06', '14:12', '10:24', '20:08', '18:46']

    assert date == reference_date, 'The date is wrong!'
    assert time == reference_time, 'The time is wrong!'

    
def test_number():
    date, time, animal, number = animals.read_animals('animals.txt')
    
    ref_numbers = [36, 25, 26, 31, 20]
    
    assert number == ref_numbers, 'The numbers are wrong!'
    
    
def test_mean1():
    l = [1]
    assert animals.mean(l) == 1

   
def test_mean2():
    l = [3, 4, 5]
    assert animals.mean(l) == 4

    
def test_mean3():
    l = []
    assert animals.mean(l) == 0

    
def test_mean4():
    l = [3.0, 4.0, 5.0]
    assert animals.mean(l) == 4

    
def test_mean5():
    l = [2, 3, 5, 4]
    assert animals.mean(l) == 3.5

    
def test_mean_animals_sighted1():
     f = 'animals.txt'
     a = 'Grizzly'
     
     m = animals.mean_animals_sighted(f, a)
     
     assert m == 36

     
def test_mean_animals_sighted2():
     f = 'animals.txt'
     a = 'Elk'
     
     m = animals.mean_animals_sighted(f, a)
     
     assert m == 25.5


import nose
@nose.tools.raises(ValueError)
def test_mean_animals_sighted3():
     f = 'animals.txt'
     a = 'Pangolin'
     
     m = animals.mean_animals_sighted(f, a)
     