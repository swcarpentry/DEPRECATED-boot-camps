
def sum_list(numbers):
    """ Sum a list of numbers and return the total.

    Arguments:
    numbers -- a list of numbers.

    """
    total = 0
    for num in numbers:
        total = total + num
    return total


numbers = [1, 2, 3, 4]
def test_sum_list(numbers):
    """ Test sum_list.

    Arguments:
    numbers -- a list of numbers.
    """
    assert sum(numbers) == 10, "mean of [1, 2, 3, 4, 5] is not 10"

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


numbers = [1, 2, 3, 4]
def test_calc_mean(numbers):
    """ Test calc_mean.

    Arguments:
    numbers -- a list of numbers.
    """
    assert calc_mean(numbers) == 2.5, "mean of [1, 2, 3, 4] is not 2.5"


test_calc_mean(numbers)
