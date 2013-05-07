from dnautils import antiparallel
from nose.tools import assert_raises

def test_antiparallel_a():
    assert antiparallel('A') == 'T'
    
def test_antiparallel_c():
    assert antiparallel('C') == 'G'

def test_antiparallel_t():
    assert antiparallel('T') == 'A'

def test_antiparallel_g():
    assert antiparallel('G') == 'C'

def test_antiparallel_gggg():
    assert antiparallel('GGGG') == 'CCCC'

def test_antiparallel_gtca():
    assert antiparallel('GTCA') == 'TGAC'

def test_antiparallel_empty_string():
    assert antiparallel('') == ''

def test_123():
    assert_raises(ValueError, antiparallel, 123)
