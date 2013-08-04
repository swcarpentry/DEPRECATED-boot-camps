def divide_it(x, y):
    try:
        out = x / y
    except:
        print '   Divide by zero!'
        out = None
    return out


divide_it(4,2)
divide_it(4,0)

