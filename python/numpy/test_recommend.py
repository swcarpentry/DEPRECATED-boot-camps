import nose.tools as nt

import recommend
from data import raw_scores

def test_load_data1():
    ratings = recommend.extract_ratings(raw_scores)
    nt.assert_equal(ratings.ndim, 2)

def test_load_data2():
    ratings = recommend.extract_ratings(raw_scores)
    nt.assert_equal(ratings.shape[0], len(raw_scores))

def test_load_data3():
    ratings = recommend.extract_ratings(raw_scores)
    nt.assert_equal(ratings.shape[1], 6)

def test_load_data4():
    ratings = recommend.extract_ratings(raw_scores)
    nt.assert_equal(ratings.sum(), 106.5)
