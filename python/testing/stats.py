def mean(numlist):
    """Calculate the arithmetic mean of a list of numbers in numlist"""
    try:
        total = sum(numlist)
        length = len(numlist)
    except TypeError:
        raise TypeError("The number list was not a list of numbers.")
    except:
        print "There was a problem evaluating the number list."
    return float(total)/length

def std(vals):
    """Computes the standard deviation from a list of values."""
    pass
