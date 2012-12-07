
pdfs/%.2d.pdf : data/%.2d.dat pdfs
	sh plot-2d.sh $< $@

pdfs/%.3d.agreement.pdf : data/%.3d.agreement.dat pdfs
	sh plot-3d-agreement.sh $< $@

pdfs/%.3d.expertise.pdf : data/%.3d.expertise.dat pdfs
	sh plot-3d-expertise.sh $< $@

pdfs/%.3d.interest.pdf : data/%.3d.interest.dat pdfs
	sh plot-3d-interest.sh $< $@

# Use wildcard function to get all .dat files and patsubst to
# convert this to a list of .pdf files.
all_dat_files := $(patsubst data/%.dat, pdfs/%.pdf, $(wildcard data/*.dat))

# Use the value of the above function to build all the pdfs.
all : $(all_dat_files)

pdfs :
	mkdir pdfs

clean :
	rm -rf pdfs/*
