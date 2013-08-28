def mean(vals):
    """Computes the mean from a list of values."""
    try:
        total = float(sum(vals))
        length = len(vals)
    except TypeError:
        raise TypeError("The list was not numbers.")
    except:
        print "Something unknown happened with the list."
    return total/length

def median(vals):
    """Computes the median from a list of values."""
    pass

def mode(vals):
    """Computes the mode from a list of values."""
    pass

def std(vals):
    """Computes the standard deviation from a list of values."""
    n = len(vals)
    if n == 0:
        return 0.0
    mu = sum(vals) / n
    if mu == 1e500:
        return NotImplemented
    var = 0.0
    for val in vals:
        var = var + (val - mu)**2
    return (var / n)**0.5

def var(vals):
    """Computes the variance from a list of values."""
    pass
