# Solutions to exercises in Write Code for People I

* * * *
##Intro to Variables

**Short Exercises**

Can you add the extra space between my last and first name?

***Solution***
```
In [3]: full_name = first_name + ' ' + last_name
In [4]: print(full_name)
Out[4]: 'PaulWilson'
```


* * * *
## Lists

**Short Exercise**

* What is the difference between the following?
   * `voltage_list[0:0]`
   * `voltage_list[0:1]`
   * `voltage_list[:1]`

* What does `current_list[-3:]` give?

* What does `voltage_list[::2]` mean?

**Solution**

* What is the difference between the following?
   * returns an empty list 
   * returns a list containing the first value
   * returns a list containing the first value
   
```
In [7]: voltage_list[0:0]
Out[7]: []

In [8]: voltage_list[0:1]
Out[8]: [-2.0]

In [9]: voltage_list[:1]
Out[9]: [-2.0]
```
   
* What does `current_list[-3:]` give?

```
In [14]: voltage_list[-3:]
Out[14]: [0.0, 1.0, 2.0]
```

* What does `voltage_list[::2]` mean?
Starts at the begining of the list and goes to the end by advancing two 
steps each time. Returns a list containing these values
```
In [15]: voltage_list[::2]
Out[15]: [-2.0, 0.0, 2.0]
```

* * * *
## Append and Extend

**Short Exercises**

Make a new list:

```
In [16]: exercise_list = [1,1,2,3,5,8]
```

What is the difference between the following?
 * `exercise_list.append([13,21,34])`
 * `exercise_list.extend([13,21,34])`

***Solution***

* `append` adds the list [13,21,34] onto the end of the list
```
Out[18]: [1, 1, 2, 3, 5, 8, [13, 21, 34]]
```

* `extend` adds the values 13, 21, 34 onto the end of the list
```
Out[23]: [1, 1, 2, 3, 5, 8, 13, 21, 34]
```


* * * *
##Tuples

**Short Exercise**

Display the second element of a tuple with two different slices.

***Solution***

```
In [1]: tup = ("red", "white", "blue")

In [2]: tup[1:2]
Out[2]: ('white',)

In [3]: tup[1::2]
Out[3]: ('white',)
```

* * * * 