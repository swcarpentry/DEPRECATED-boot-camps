def scan(lst):
    total = 0
    last = None
    if len(lst) == 0:
        return total
    for item in lst:
        if last is None or item > last:
            total += item
            last = item
    return total
