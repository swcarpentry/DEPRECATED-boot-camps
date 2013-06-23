def divide_it(x, y):
    try:
        out = x / y
    except ZeroDivisionError:
        print '   Divide by zero!'
        out = None
    except TypeError:
        print ' Provide a number!'
        out = None    
    return out

def divide_it1(x,y):
    return x/y

divide_it(4,2)
divide_it(4,0)
divide_it('a',2)
