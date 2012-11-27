#This Makefile runs test_functions.py if calc_mean.py has been updated.
#It will also run calc_mean.py on student_grades1.dat if student_grades1.dat has been updated

test_functions.py : calc_mean.py
	nosetests test_functions.py

calc_mean.py : student_grades1.dat
	python calc_mean.py student_grades1.dat


