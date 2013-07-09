
# Good Programming Practice - Documenting Code

In the [last session](2_functions_and_modules.md) you learned how to package code into functions and to package functions into modules (also called libraries). Functions and modules let you easily design, write and package your code so that it is easy to understand and easily reusable. However, to share the code, and really understand what it works, you need to add documentation.

You have already seen documentation using python "help()". For example, lets look at the documentation for the "string" module that we used in the [first session](1_lists_and_dictionaries.md).

    $ ipython
    $ import string
    $ help(string)
    Help on module string:
    
    NAME
        string - A collection of string operations (most are no longer used).
    
    FILE
        /System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/string.py
    
    MODULE DOCS
        http://docs.python.org/library/string
    
    DESCRIPTION
        Warning: most of the code you see here isn't normally used nowadays.
        Beginning with Python 1.6, many of these functions are implemented as
        methods on the standard string object. They used to be implemented by
        a built-in module called strop, but strop is now obsolete itself.

Lets compare this to the documentation for the "checkmain.py" script that we saw in the last session.

    $ import checkmain
    $ help(checkmain)
    Help on module checkmain:
    
    NAME
        checkmain
    
    FILE
        /Users/chris/Work/Teaching/swcarpentry/course/checkmain.pyc
    
    FUNCTIONS
        addArrays(x, y)

Not great... It is very important when programming in any language that we provide full documentation for all of the functions and modules. In python, this is achieved by adding documentation strings to each part of the script. These are strings that are placed at the beginning of the function or module.

    $ def documentedFunction(a):
    $     """Here is the documentation string for this function"""
    $     return a
    $
    $ help(documentedFunction)
    
    Help on function documentedFunction in module __main__:
    
    documentedFunction(a)
        Here is the documentation string for this function

We can do the same thing for the [checkmain.py](checkmain.py) script;

    """checkmain is a simple python script to demonstrate
       hiding the code if the script is imported as a module"""

    def addArrays(x, y):
        """This function adds together each element of the two
           passed lists, returning the result in the returned list."""
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

We now get better documentation when using help()

    $ ipython
    $ import checkmain
    $ help(checkmain)
    Help on module checkmain:
    
    NAME
        checkmain
    
    FILE
        /Users/chris/Work/Teaching/swcarpentry/course/checkmain.py
    
    DESCRIPTION
        checkmain is a simple python script to demonstrate
        hiding the code if the script is imported as a module
    
    FUNCTIONS
        addArrays(x, y)
            This function adds together each element of the two
            passed lists, returning the result in the returned list.

## Exercise

### Exercise 3

Edit your [morse.py](2b/example/morse.py) script and add documentation strings for the module and also for all of the functions.

If you are really stuck then there is a completed example script in [3/example/morse.py](3/example/morse.py)

Make sure that you commit your edited script to your Git repository.

    $ git pull
    $ git commit -a
    $ git push

# [Previous](2_functions_and_modules.md) [Up](python_and_good_programming_practice.md) [Next](4_object_orientation.md) 
