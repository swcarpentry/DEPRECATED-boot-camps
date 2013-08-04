## Test-driven development


[Test driven development](http://www.amazon.com/Test-Driven-Development-By-Example/dp/0321146530) (TDD), proposed by Kent Beck, is a philosophy and approach where we write code by *writing the tests first*, then write the code to make the tests pass. If a new feature is needed, another test is written and the code is expanded to meet this new use case. This continues until the code does what is needed. This can be summarised as red-green-refactor:

 * Red - write tests based on requirements. They fail as there is no code!
 * Green - write/modify code to get tests to pass.
 * Refactor code - clean it up.
 
The above steps are usually repeated in a number of iterations as you gather and refine the requirements.

By writing tests first, we're forced to think about what our code should do. In contrast, in writing our code then tests, we risk testing what the code actually *does*, rather than what it *should* do.

TDD operates on the YAGNI principle (You Ain't Gonna Need It) to avoid developing code for which there is no need.

We'll do some simple TDD using an example with different animals and their behaviours. Let's start with some tests (saving them in [animals_0.py](python-code/animals/animals_0.py)):

    def test_moves():
        assert Animal('owl').move() == 'fly'
        assert Animal('cat').move() == 'walk'
        assert Animal('fish').move() == 'swim'

    def test_speaks():
        assert Animal('owl').speak() == 'hoot'
        assert Animal('cat').speak() == 'meow'
        assert Animal('fish').speak() == ''

If we run `nosetests animals_0.py` we'll get an error - we need to add an 'Animal' class that meets our requirements. Just in order to be able to keep track of how our TDD develops, let's save the code in a different file [animals_1.py](python-code/animals/animals_1.py).

    class Animal:
        animal_defs = {'owl':{'move':'fly','speak':'hoot'},
                   'cat':{'move':'walk','speak':'meow'},
                   'fish':{'move':'swim','speak':''}}

        def __init__(self, name):
            self.name = name

        def move(self):
            return self.animal_defs[self.name]['move']
        
        def speak(self):
            return self.animal_defs[self.name]['speak']


    def test_moves():
        assert Animal('owl').move() == 'fly'
        assert Animal('cat').move() == 'walk'
        assert Animal('fish').move() == 'swim'

    def test_speaks():
        assert Animal('owl').speak() == 'hoot'
        assert Animal('cat').speak() == 'meow'
        assert Animal('fish').speak() == ''


Our tests helped us in writing our class - we knew what we wanted in an animal definition. But now we want our animals to become more active ([animals_2.py](python-code/animals/animals_2.py)):

    from random import random

    class Animal:
        animal_defs = {'owl':{'move':'fly','speak':'hoot'},
                   'cat':{'move':'walk','speak':'meow'},
                   'fish':{'move':'swim','speak':''}}
                   
        def __init__(self, name):
            self.name = name

        def move(self):
            return self.animal_defs[self.name]['move']
        
        def speak(self):
            return self.animal_defs[self.name]['speak']


    def test_moves():
        assert Animal('owl').move() == 'fly'
        assert Animal('cat').move() == 'walk'
        assert Animal('fish').move() == 'swim'

    def test_speaks():
        assert Animal('owl').speak() == 'hoot'
        assert Animal('cat').speak() == 'meow'
        assert Animal('fish').speak() == ''

    def test_dothings_list():
        """ Test that the animal does the same number of things as the number of hour-times given."""
        times = []
        for i in xrange(5):
            times.append(random() * 24.)
        for a in ['owl', 'cat', 'fish']:
            assert len(Animal(a).dothings(times)) == len(times)

    def test_dothings_with_beyond_times():
        for a in ['owl', 'cat', 'fish']:
            assert Animal(a).dothings([-1]) == ['']
            assert Animal(a).dothings([25]) == ['']

    def test_nocturnal_sleep():
         """ Test that an owl is awake at night."""
        night_hours = [0.1, 3.3, 23.9]
        noct_behaves = Animal('owl').dothings(night_hours)
        for behave in noct_behaves:
            assert behave != 'sleep'

Now we have to write the function (`dothings`) for which we just implemented the tests ([animals_3.py](python-code/animals/animals_3.py)).

    from random import random

    class Animal:
        animal_defs = {'owl':{'move':'fly','speak':'hoot'},
                   'cat':{'move':'walk','speak':'meow'},
                   'fish':{'move':'swim','speak':''}}
                   
        def __init__(self, name):
            self.name = name

        def move(self):
            return self.animal_defs[self.name]['move']
        
        def speak(self):
            return self.animal_defs[self.name]['speak']

    def dothings(self, times):
        """ A method which takes a list of times (hours between 0 and 24) and returns a list of what the animal is (randomly) doing. - Beyond hours 0 to 24: the animal does: ""  """
        out_behaves = []
        for t in times:
            if (t < 0) or (t > 24):
                out_behaves.append('')
            elif ((self.name == 'owl') and (t > 6.0) and (t < 20.00)):
                out_behaves.append('sleep')
            else:
                out_behaves.append(self.animal_defs[self.name]['move'])
        return out_behaves


    def test_dothings_list():
        """ Test that the animal does the same number of things as the number of hour-times given."""
        times = []
        for i in xrange(5):
            times.append(random() * 24.)
        for a in ['owl', 'cat', 'fish']:
            assert len(Animal(a).dothings(times)) == len(times)

    def test_dothings_with_beyond_times():
        for a in ['owl', 'cat', 'fish']:
            assert Animal(a).dothings([-1]) == ['']
            assert Animal(a).dothings([25]) == ['']

    def test_nocturnal_sleep():
         """ Test that an owl is awake at night."""
        night_hours = [0.1, 3.3, 23.9]
        noct_behaves = Animal('owl').dothings(night_hours)
        for behave in noct_behaves:
            assert behave != 'sleep'

    if __name__ == '__main__':
       
        c = Animal('cat')
        o = Animal('owl')
        f = Animal('fish')

        times = []
        for i in xrange(10):
            times.append(random() * 24.)
        times.sort()
    
        c_do = c.dothings(times)
        o_do = o.dothings(times)
        f_do = f.dothings(times)
        
        #Let's see what the animals are actualy doing 
         
        for i in xrange(len(times)):
            print "time=%3.3f cat=%s owl=%s fish=%s" % (times[i], c_do[i], o_do[i], f_do[i])



Previous: [Testing in practice](RealWorld.md) Next: [Conclusions and further information](Conclusion.md)
