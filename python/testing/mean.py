def mean(vals):
    try :
        total = sum(vals)
        length = len(vals)
    except TypeError :
        raise TypeError("The list contained non-numeric elements.")
    except :
        print "Something unknown happened with the list."
    return total/length
