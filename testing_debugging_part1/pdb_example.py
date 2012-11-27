import sys
'''
This code is intentionally uncommented. You should never write code that looks like this!
'''


def function1(a, b, c):
    x = a
    y = b
    z = c
    return x, y, z

def function2(n):
    x = range(n)
    y = 2*x
    z = x
    return x, y, z

if __name__ == "__main__":
    a = float(sys.argv[1])
    b = float(sys.argv[2])
    c = float(sys.argv[3])

    x, y, z = function1(a, b, c)
    
    n = int(sys.argv[4])

    x, y, z= function2(n)
    z[n - 1] = n+4

    if x == z:
        print 'x = z'
