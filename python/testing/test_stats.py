from stats import mean

def test_mean():
    """Test some standard behavior of the mean() function."""
    assert(mean([2, 4]) == 3)
    assert(mean([2,4,6]) == 4)
    assert(mean([2, 4, 6, 8]) == 5)

def test_negative_mean():
    """Test standard behavior of the mean() function with negative numbers."""
    assert(mean([-4, -2]) == -3)
    assert(mean([-2, 1, 4]) == 1)

def test_float_mean():
    """Test standard behavior of the mean() function with floats numbers."""
    assert(mean([1, 2]) == 1.5)
    assert(mean([2.0, 4.0, 6.0]) == 4.0)
    assert(mean([2.5, 4.5, 6.5]) == 4.5)

def test_string_mean():
    try:
        mean(['hello','world!'])
    except TypeError:
        pass

test_mean()
test_negative_mean()
test_float_mean()
test_string_mean()
