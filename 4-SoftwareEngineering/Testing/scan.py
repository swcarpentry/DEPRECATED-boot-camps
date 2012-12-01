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
