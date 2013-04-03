from cdna import complement

def test_complement_a():
    assert complement('A') == 'T'
    
def test_complement_c():
    assert complement('C') == 'G'

def test_complement_t():
    assert complement('T') == 'A'

def test_complement_g():
    assert complement('G') == 'C'

def test_complement_gatc():
    assert complement('GATC') == 'CTAG'

def test_complement_empty_string():
    assert complement('') == ''
