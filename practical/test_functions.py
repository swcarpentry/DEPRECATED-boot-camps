import calc_mean
from calc_mean import get_mean, read_file

def test_mean_functionality():
    x = [1, 2, 3.0, 5, 1.0, 0]
    assert get_mean(x) == 2.0, 'Failed mean_test_functionality, mean([1, 2, 3.0, 5. 1.0, 0]) did not equal 2.0'

def test_mean_negative():
    x = [-1, -2, 4, -5, -9]
    assert get_mean(x) == -2.6, 'Failed mean_test_negative, mean([-1, -2, 4, -5, -9]) did not equal -2.6'

def test_mean_invalid_list_entry():
    x = [1, 2, 3, '4', 5]
    mean = get_mean(x)
    assert mean == 'Error', 'Failed mean_test_invalid_list_entry, program failled to exit with an assertion error when a string was encountered'

def test_read_file_functionality():
    grades = read_file('student_grades_test1.dat')
    assert grades == [91.0, 90.0, 92.0], 'Failed read_file_test_functionality, grades did not equal [91.0, 90.0, 92.0]'

def test_read_file_blank_line():
    grades = read_file('student_grades_blank.dat')
    assert grades == 'Error', 'Failed read_file_test_blank_line, program failed to exit with an assertion error when a blank line was encountered'
