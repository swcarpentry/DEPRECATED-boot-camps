# Use POSIX-compatible processing [target-rules,usage]
.POSIX:

# Define macros (variables) for more general rules [macros]
CSV_SOURCE = kepler_planets.csv exoplanet.eu_catalog.csv
PYTHON = python
HTML = exoplanet_figures.html

# Declare suffixes for inference rules [target-rules]
.SUFFIXES: .csv .png

# Don't clean up prerequisites [target-rules]
.PRECIOUS:

# The first target is the default [operands]
all: $(HTML)

# Inference rule building %.png from %.csv [inference-rules]
.csv.png:
	$(PYTHON) make_plots.py $^

# Build HTML output based on PNG input
$(HTML): $(CSV_SOURCE:.csv=.png)
	$(PYTHON) make_html.py $@ $^

# Cleanup any automatically-generated files
clean:
	rm -f *.png $(HTML)

# References:
# [operands]:
#   http://pubs.opengroup.org/onlinepubs/9699919799/utilities/make.html#tag_20_76_05
# [target-rules]:
#   http://pubs.opengroup.org/onlinepubs/9699919799/utilities/make.html#tag_20_76_13_04
# [macros]:
#   http://pubs.opengroup.org/onlinepubs/9699919799/utilities/make.html#tag_20_76_13_05
# [inference-rules]:
#   http://pubs.opengroup.org/onlinepubs/9699919799/utilities/make.html#tag_20_76_13_06
# [usage]:
#   http://pubs.opengroup.org/onlinepubs/9699919799/utilities/make.html#tag_20_76_16
