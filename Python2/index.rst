Python Data Structures
======================

**Updated and presented by: Tracy Teal**

Adapted from Software Carpentry materials
http://software-carpentry.org/

Starting an iPython notebook
Mac
~/anaconda/bin/ipython notebook --pylab=inline
PC
./run-in-vm.sh



Python Lists
------------
SWC Tutorial: http://software-carpentry.org/4_0/python/lists/


Collections let us store many values together

The most common way we do this is with a list

We create a list in Python with  ::

   listname = ['a', 'b', 'c', 'd']

example ::

   gases = ['He', 'Ne', 'Ar', 'Kr']

The index for lists starts with 0 instead of 1, so the first item in a list
is item 0 ::

   gases[0] would return 'He'
   gases[1] would return 'Ne'

You can also get items from the end of the list ::

   gases[-1] would return 'Kr'
   gases[-2] would return 'Ar'

Use  ::

   len(listname) 

to get the length of the list, or how many values are in the list.

You can also change a list after you make it.  If you want to change, 
say Kr to K you can do ::

   gases[3] = 'K'

Now your list will be:  ['He', 'Ne', 'Ar', 'K']

Lists can store values of many kinds, even other lists ::

   helium = ['He', 2]
   neon = ['Ne', 8]
   gases = [helium, neon]

Now if you want do something to every item in the list, you can use a loop ::

   gases = ['He', 'Ne', 'Ar', 'Kr']
   i = 0
   while i < len(gases):
      print gases[i]
      i += 1

This will print out each of the gases.

A better way to do this would be to use a 'for' loop ::

    for i in gases:
       print i

This will also print out each gas.

You can use if statements to see if something in the list is true ::

   if 'Pu' in gases:
      print 'But plutonium is not a gas'
   else:
      print 'The universe is well ordered'



Some useful list methods
-------------------------

You can append items to the list::

   gases.append('H')

You can print out how many of something is in the list::
  
   print gases.count('He')

You can print where the item is in the list::
 
   print gases.index('Ne')



Dictionaries
------------
SWC Tutorial: http://software-carpentry.org/4_0/setdict/dict/


Dictionaries have key, value pairs.  Here is an example of a dictionary.

>>> tel = {'jack': 4098, 'sape': 4139}
>>> tel['guido'] = 4127
>>> tel
{'sape': 4139, 'guido': 4127, 'jack': 4098}
>>> tel['jack']
4098
>>> del tel['sape']
>>> tel['irv'] = 4127
>>> tel
{'guido': 4127, 'irv': 4127, 'jack': 4098}
>>> tel.keys()
['guido', 'irv', 'jack']
>>> 'guido' in tel
True

You can loop over items in a dictionary the same way you can over items in a 
list.  ::

    for keys in tel:
      print tel[keys]

If you want to mix some text in with your printing ::

    for keys in tel:
      print 'This is the number', tel[keys]


