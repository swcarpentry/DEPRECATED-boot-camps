
def addArrays(x, y):
    z = []
    for i in range(0,len(x)):
       z.append( x[i] + y[i] )
    
    return z
    
if __name__ == "__main__":
    # Don't run this code if this script is being
    # imported as a module 

    a = [ 1, 2, 3, 4 ]
    b = [ 5, 6, 7, 8 ]

    c = addArrays(a, b)
    print( c )

