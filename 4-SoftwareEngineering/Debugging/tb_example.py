def check_n(n):
    assert 0 <= n
    assert isinstance(n, int)


def fib(n):
    if 0 <= n <= 1:
        return n
    elif 1 < n:
        return fib(n-1) + fib(n-2)
    else:
        check_n(n)


def main():
    a = 4.2
    fib(a)


if __name__ == "__main__":
    main()
