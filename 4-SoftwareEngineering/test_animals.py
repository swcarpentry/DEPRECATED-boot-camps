import animals

def test_read_animals_date():
    reference_date = ['2011-04-22', '2011-04-23', '2011-04-23', '2011-04-23', 
                      '2011-04-23']

    date, time, animal, number = animals.read_animals('animals.txt')

    assert date == reference_date, 'Date is wrong'

def test_read_animals_time():
    reference_time = ['22:06', '14:12', '10:24', '20:08', '18:46']

    date, time, animal, number = animals.read_animals('animals.txt')

    assert time == reference_time, 'Time is wrong'
