#tex stuff
TEX = $(wildcard *.tex)
PDF = $(UI:.tex=.pdf)
SLAG = $(wildcard *.out *.log *.aux *.nav *.snm *.toc) $(HTMLTARGET)

#landslide stuff
RST=$(wildcard *.rst)
HTMLTARGET = $(RST:.rst=.html)

#hack to cope with trailing spaces on src files
TRAILING=$(TEX:.tex=.tex.trailing) $(RST:.rst=.rst.trailing)

DPI=600
SVG=$(wildcard *.svg)
PNGS=$(SVG:.svg=-$(DPI)dpi.png)

all: $(PNGS)
	@echo See $^
.PHONY: all

svgpng=inkscape \
	--without-gui --export-area-page --export-dpi=$1 --file=$2 --export-png=$3

flatten=convert -flatten +matte $1 $2

.PHONY: clean
.IGNORE: clean

help:
	@echo By default, builds PNGs at 600 dpi by default.
	@echo For instance, call "make DPI=200"

all: $(PDF) $(HTMLTARGET)

clean-build:
	-rm $(SLAG) $(TRAILING)

clean: clean-build

trailing-spaces: $(TRAILING)

show: all
	see $(HTMLTARGET)

show-chromium: all
	chromium $(HTMLTARGET)

%.html: %.rst
	landslide -t ./themes -i $< -d $@

$(PDF): $(TEX)

%.pdf: %.tex
	pdflatex $<


%.trailing: %
	sed -i 's/[ \t]*$$//' $<
	touch $@ # stamp to avoid re seding
	@echo Removal of trailing spaces: should be taken care of by vim.
