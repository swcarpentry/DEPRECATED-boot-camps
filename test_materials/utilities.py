
def sum_list(numbers):
    """ Sum a list of numbers and return the total.

    Arguments:
    numbers -- a list of numbers.

    """
    total = 0
    for num in numbers:
        total = total + num
    return total

def test_sum_list(numbers):
    """ Test sum_list.

    Arguments:
    numbers -- a list of numbers.

    """
    assert sum(numbers) == 10, "sum of %s is not 10" % str(numbers)
numbers = [1, 2, 3, 4]
test_sum_list(numbers)

def calc_mean(numbers):
    """ Calculate the mean of a list of numbers.

    Arguments:
    numbers -- a list of numbers.

    """
    total = 0
    for num in numbers:
        total = total + num
    return(total / float(len(numbers)))

def test_calc_mean(numbers):
    """ Test calc_mean.

    Arguments:
    numbers -- a list of numbers.

    """
    assert calc_mean(numbers) == 2.5, "mean of %s is not 2.5" % str(numbers)
numbers = [1, 2, 3, 4]
test_calc_mean(numbers)
