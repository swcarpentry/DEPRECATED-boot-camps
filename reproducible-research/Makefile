#tex stuff
TEX = $(wildcard *.tex)
PDF = $(UI:.tex=.pdf)
SLAG = $(wildcard *.out *.log *.aux *.nav *.snm *.toc) $(HTMLTARGET)

#landslide stuff
RST=$(wildcard *.rst)
HTMLTARGET = $(RST:.rst=.html)

#hack to cope with trailing spaces on src files
TRAILING=$(TEX:.tex=.tex.trailing) $(RST:.rst=.rst.trailing)

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
