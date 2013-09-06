
# Program Design as Exemplified in NumPy

In the last session you saw how to use classes and objects to write software that can be used safely and easily by other developers. It is good practice to write software that can be used easily by other people. The best way to learn how to do this is to read and try to use other people's software. You will learn to improve your code by working out what they did right, and also by getting annoyed by the things that they did wrong.

NumPy is a complete, object-orientated module for fast maths (numerics) in Python. All of the code needed to manipulate matrices and vectors, calculate random numbers and perform polynomial fitting is obviously very complex. This complexity is hidden behind a set of classes / objects that are provided by NumPy to provide a simple, easy-to-use interface for users such as us.

First, load the NumPy module.

    $ ipython
    $ import numpy

The first thing to do when looking at some new code is to try to read the documentation...

    $ help(numpy)
    Help on package numpy:
    
    NAME
        numpy
    
    FILE
        /System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python/numpy/__init__.py
    
    DESCRIPTION
        NumPy
        =====

The NumPy documentation is very complete and useful. You will also see that it points you towards a reference guide that is held at [http://scipy.org](http://scipy.org).

NumPy is built from a series of submodules, that can be seen listed near the end of the help. Modules include fft (fast fourier transform), linalg (linear algebra), polynomial (polynomial fitting etc.) and math (access to standard math functions such as cosine, sine etc.).

You can get the help on the fft submodule using

    $ help(numpy.fft)

All of numpy's modules are loaded when you type 

    $ import numpy

So, for example, the "cos" function in the math module can be used via

    $ numpy.math.cos(0)
    1.0

However, you can also import a single module, e.g. to import just the math module, type

    $ from numpy import math

Now, the cos function can be used via

    $ math.cos(0)
    1.0

You can also rename modules when you load them. This can be used to create a quick shorthand, e.g. renaming numpy as "np"

    $ import numpy as np
    $ np.math.cos(0)
    1

or

    $ from numpy import math as nmath
    $ nmath.cos(0)
    1

## Exercise

### Exercise 5a

NumPy is famous for fast matrix algebra. Look at the ipython help for the NumPy matrix class;

    # help(numpy.matrix)

First, use the help to work out how to create the 2D identity matrix (1,0 ; 0,1)?
Next, use the same technique to create a matrix corresponding to the vector (1,1).

The benefit of the NumPy matrix class, is that multiplication of matrix objects will follow the rules of matrix maths. Verify this by multiplying your vector by the identity matrix. The resulting vector should be the same.

Now, create a 2D matrix with values (2,0 ; 0,2). This should scale your vector by two times. Verify that multiplying your vector by this matrix returns the vector (2,2).

Next, while you can multiply a vector by a matrix, you cannot multiply a matrix by a vector. Verify that a python exception will be raised when you multiply the identity matrix by the vector.

Finally, turn the above four tests into a nose-tests script to test numpy.matrix. Note that to test for equality of two matricies, you need to use;

    $ assert( (matrix1 == matrix2).all() )

If you are stuck, there is an example completed script in (5/example/matrix_test.py)[5/example/matrix_test.py]
