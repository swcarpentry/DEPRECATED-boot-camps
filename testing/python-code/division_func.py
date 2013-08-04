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

DIVISION_RESULTS = {'4/2': 2, '6/3':2, '10/2':5}