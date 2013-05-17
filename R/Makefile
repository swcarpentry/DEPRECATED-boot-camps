all: figs/slopes_AsiaVsAmericas.pdf

results/gCoef.txt results/gCoef.rds: data/gapminderDataFiveYear.txt code/01_countrySpecificInterceptSlope.R
	R CMD BATCH code/01_countrySpecificInterceptSlope.R

figs/slopes_AsiaVsAmericas.pdf results/slopes_AsiaVsAmericas.txt: results/gCoef.rds code/03_slopeComparisonAsiaVsAmericas.R
	R CMD BATCH code/03_slopeComparisonAsiaVsAmericas.R