from rectangle import overlap

def test_coincident():
    a = [0, 2, 0, 2]
    b = a
    expected = a
    actual = overlap(a, b)
    assert expected == actual

def test_a_encloses_b():
    a = [0, 3, 0, 3]
    b = [1, 2, 1, 2]
    expected = b
    actual = overlap(a, b)
    assert expected == actual

def test_b_encloses_a():
    a = [1, 2, 1, 2]
    b = [0, 3, 0, 3]
    expected = a
    actual = overlap(a, b)
    assert expected == actual

def test_a_top_right_b():
    a = [3, 6, 3, 6]
    b = [0, 4, 0, 4]
    expected = [3, 4, 3, 4]
    actual = overlap(a, b)
    assert expected == actual

def test_a_left_of_b():
    a = [0, 2, 0, 2]
    b = [3, 4, 3, 4]
    expected = None
    actual = overlap(a, b)
    assert expected == actual
