#This Makefile runs test_functions.py if calc_mean.py has been updated.
#It will also run calc_mean.py on student_grades1.dat if student_grades1.dat has been updated

#run nosetests always
test : test_functions.py
	nosetests test_functions.py

#Run nosetests if calc_mean.py or test_functions.py are updates
nosetests.xml : test_functions.py calc_mean.py
	nosetests test_function.py --with-xunit

# If calc_mean.py is updated then mark test_functions.py as updated.
test_functions.py : calc_mean.py
	touch test_functions.py


#Run calc_mean on student_gradesN.dat
%.dat.avg : %.dat calc_mean.py
	python calc_mean.py $<>$@


