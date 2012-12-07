
pdfs/%.pdf : data/%.2d.dat pdfs
	sh plot-2d.sh $< $@

pdfs/%.2d.pdf : data/%.2d.dat pdfs
	sh plot-2d.sh $< $@

pdfs/%.3d.agreement.pdf : data/%.3d.agreement.dat pdfs
	sh plot-3d-agreement.sh $< $@

pdfs/%.3d.expertise.pdf : data/%.3d.expertise.dat pdfs
	sh plot-3d-expertise.sh $< $@

pdfs/%.3d.interest.pdf : data/%.3d.interest.dat pdfs
	sh plot-3d-interest.sh $< $@

pdfs :
	mkdir pdfs

clean :
	rm -rf pdfs/*
