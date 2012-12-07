
pdfs/%.pdf : data/%.2d.dat pdfs
	sh plot-2d.sh $< $@

pdfs :
	mkdir pdfs

clean :
	rm -rf pdfs/*
