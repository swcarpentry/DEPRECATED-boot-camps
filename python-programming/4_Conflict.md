#Python and git - version controlling your code

We are now going to create a script which calculates the area of a square. 
This script is going to contain a function called area that calculates the area
and prints it to screen. It will also contain a __main__ if statement that will
let you test your code before you import it into other scripts.

You will be working in pairs - one of you will create a version that has a function
defined like this:

```python
def area(x1, y1, x2, y2):
```

and the other one will create a function that is defined like this:

```python
def area(point1, point2):
```
where each of the points are a list with the x and the y coordinate.

So - one of you has the coordinates as just a list of numbers, and the other one
has this a set of two lists, each list containing one coordinate set.

Note: in both cases you can assume that the values of the second coordinate pair is higher
than that of the first. 

Name your python script area.py.


## Solutions(s)

The two different versions of this script look like this:

```python
def area(x1, y1, x2, y2):
    xdiff = x2 - x1
    ydiff = y2 - y1
    area = xdiff * ydiff
    return area

if __name__ == "__main__":
    print area(2,2,4,4)

```

```python
def area(point1, point2):
    xdiff = point2[0] - point1[0]
    ydiff = point2[1] - point1[1]
    area = xdiff * ydiff
    return area

if __name__ == "__main__":
    firstpoint = [2,2]
    secondpoint = [4,4]
    print area(firstpoint, secondpoint)
```

