from division_func import divide_it
from division_func import DIVISION_RESULTS

def test_4_2():
    assert divide_it(4,2) ==  DIVISION_RESULTS['4/2']
def test_6_3():
    assert divide_it(6,3) == DIVISION_RESULTS['6/3']
def test_10_2():
    assert divide_it(10,2) == DIVISION_RESULTS['10/2']

test_4_2()
test_6_3()
test_10_2()