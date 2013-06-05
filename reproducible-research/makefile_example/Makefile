
TEX=$(wildcard *.tex)
PDF=$(TEX:.tex=.pdf)


all: $(PDF)


results/data.npy: generate_dataset.py
	mkdir -p results
	python generate_dataset.py


results/fit.npy: fit.py results/data.npy
	python fit.py


images/my_plot.png: plot.py results/fit.npy
	mkdir -p images
	python plot.py

%.pdf: %.tex images/my_plot.png
	pdflatex $<

clean:
	rm -rf images
	rm -rf results
