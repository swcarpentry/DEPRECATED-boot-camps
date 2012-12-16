all: exoplanet_figures.html

*.png: *.csv
	make_plots.py $^

exoplanet_figures.html: *.png
	make_html.py $@ $^

clean:
	rm *.png exoplanet_figures.html
