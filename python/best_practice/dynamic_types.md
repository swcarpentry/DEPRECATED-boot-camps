[Up To Schedule](../../README.md) - Return To [Write Code for People](Readme.md)

- - - -


## Types and Dynamic Typing

Like in most programming languages, things in python are *typed* — the *type* refers to the type of data and what you can do with it. For example, you can do different things with strings and numbers. Numbers can have decimal components or not, and so on. You can see the type of a variable with the `type` command.

```
In [5]: type?
Type:       type
String Form:<type 'type'>
Namespace:  Python builtin
Docstring:
type(object) -> the object's type
type(name, bases, dict) -> a new type

In [6]: type(full_name)
Out[6]: str
In [7]: type(10)
Out[7]: int

```


Python is what is known as a *dynamically typed* language. Dynamic typing means that you don't have to declare the type of a variable when you define it; python just figures it out based on how you are setting the variable. This is in contrast to *statically typed* languages, where you must say up front that a variable is going to be used for strings or numbers or whatever. There are good and bad points to both approaches.

But back to Python. Let's say you set a variable. Sometime later you can just change the type of data assigned to a variable and python is perfectly happy about that. Since it won't be obvious until (possibly much) later why that's important, I'll let you marinate on that idea for a second.

Here's an example of dynamic typing:

```
In [8]: voltage = 2
In [9]: print(voltage)
2

In [10]: type(voltage)
Out[10]: int
```

It's an `int`, which is an integer — a number with no decimal component. Let's assign a value of 2.7 (which has a decimal part) to voltage. What happens to the type?

```
In [11]: voltage = 2.7

In [12]: type(voltage)
Out[12]: float
```

Neat! That's a `float`, which does have a decimal part. You can assign a string to the variable voltage in the same way:

```python
In [13]: voltage = "2.7 volts"

In [14]: type(voltage)
Out[14]: str
```

I'll let you ruminate on the pros and cons of this construction while I change the value of voltage back to an int:

```
In [15]: voltage = 2
```

Choosing an appropriate variable type is not just a practical concern; it can also have an effect on code readability. Is this number used for calculations or only in print statements?


## On Being Precise With floats and ints

Again, the following may seem esoteric and pedantic, but it is very important. So bear with me.

Lets say you had some voltage data that looks like the following

```
0
0.5
1
1.5
2
```

If you just assigned this data individually to a variable, you'd end up with the following types:

```
0   -> int
0.5 -> float
1   -> int
1.5 -> float
2   -> int
```

But what if you wanted all of that data to be floats on its way in? You could assign the variable and then coerce it to type float:

```
In [28]: voltage = float(1)
```

But that's ugly. If you want whats otherwise an integer to be a float, add `.0` to the end:

```python
In [29]: voltage = 1.0

In [30]: type(voltage)
Out[30]: float
```

This point becomes important when we start operating on data in the next section.

## Data Operations

In this section all of the discussion in the previous section becomes
important. I don't know if I'd call this stuff fundamental to the language,
but it's pretty important and it will zing you if you aren't careful. The
takeaway is that you need to be precise with what you are doing. Let's say you
want to add some integers.

```
In [31]: a = 1

In [32]: b = 2

In [33]: c = a + b

In [34]: c
Out[34]: 3

In [38]: type(a), type(b), type(c)
Out[38]: (int, int, int)
```

So we got a value of three for the sum, which also happens to be an
integer. Any operation between two integers is another integer. Makes sense.

So what about the case where a is an integer and b is a float?

```
In [39]: a = 1

In [40]: b = 2.0

In [41]: c = a + b

In [42]: c
Out[42]: 3.0

In [43]: type(a), type(b), type(c)
Out[43]: (int, float, float)
```

You can do multiplication on numbers as well.

```
In [44]: a = 2

In [45]: b = 3

In [46]: c = a * b

In [47]: c
Out[47]: 6

In [48]: type(a), type(b), type(c)
Out[48]: (int, int, int)
```

Also division.

```
In [49]: a = 1

In [50]: b = 2

In [51]: c = a / b

In [52]: c
Out[52]: 0
```

**ZING!**

Here's why type is important. Dividing two integers returns an integer: this operation calculates the quotient and floors the result to get the answer.

If everything was a float, the division is what you would expect.

```
In [53]: a = 1.0

In [54]: b = 2.0

In [55]: c = a / b

In [56]: c
Out[56]: 0.5

In [57]: type(a), type(b), type(c)
Out[57]: (float, float, float)
```

- - - -

[Up To Schedule](../../README.md) - Return To [Write Code for People](Readme.md)

