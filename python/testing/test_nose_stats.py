from stats import mean
from nose.tools import *

def test_mean():
    """Test some standard behavior of the mean() function."""
    assert_equal(mean([2, 4]),        3)
    assert_equal(mean([2,4,6]),       4)
    assert_equal(mean([2, 4, 6, 8]),  5)

def test_negative_mean():
    """Test standard behavior of the mean() function with negative numbers."""
    assert_equal(mean([-4, -2]),    -3)
    assert_equal(mean([-2, 1, 4]),   1)

def test_float_mean():
    """Test standard behavior of the mean() function with floats numbers."""
    assert_equal(mean([1, 2]),1.5)
    assert_equal(mean([2.0, 4.0, 6.0]),4.0)
    assert_equal(mean([2.5, 4.5, 6.5]),4.5)
    assert_equal(mean([2.5, 4.5, 6.0]),4.3333)

def test_string_mean():
    assert_raises(TypeError,mean,['hello','world'])

