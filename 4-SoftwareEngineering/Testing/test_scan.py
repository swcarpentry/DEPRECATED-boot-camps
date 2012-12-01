def scan(lst):
    if len(lst) == 0:
        return 0
    total = 0
    last = None
    for item in lst:
        if last is None or last < item:
            total += item
            last = item
        elif last > item:
            total = 0
            last = None
        else:
            pass
    return total

assert scan([])                 == 0
assert scan([1, 2, 3, 4])       == 10
assert scan([2, 3, 4, 5])       == 14
assert scan([1, 1, 1, 1])       == 1
assert scan([3, 2, 5, 4])       == 0
assert scan([4, 3, 2, 1])       == 0
assert scan([3, 4, 5, 1])       == 0
assert scan([1, 2, 3, 5, 4])    == 0
assert scan([1, 2, 3, 3])       == 6
assert scan([1, 2, 4, 4])       == 7
assert scan([1, 2, 4, 4, 5])    == 12
assert scan([-3, -2, -1])       == -6
assert scan([-1, -2, -3])       == -3
assert scan([-1, -2, -3, -2, -1]) == -6
assert scan([3, 2, 1, 2, 3])      == 6
