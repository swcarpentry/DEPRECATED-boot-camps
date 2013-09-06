
import numpy as np

identity = np.matrix( [ [1,0], [0,1] ] )
scale = np.matrix( [ [2,0], [0,2] ] )
vector = np.matrix( [ [1],[1] ] )

def test_identity():
    assert( (identity*vector == vector).all() )

def test_scale():
    expected_result = np.matrix( [ [2],[2] ] )
    assert( (scale*vector == expected_result).all() )

def test_correct_maths():
    # should not raise an exception
    identity * vector

    # should raise an exception
    try:
        vector * identity
        assert(False)
    except:
        pass
