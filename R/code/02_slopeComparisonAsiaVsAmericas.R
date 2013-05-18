## use File --> Compile Notebook on me !
## or click notebook button in top right of editor pane
## will require knitr package to work
## to install: install.packages("knitr")

library(lattice)
## this path is affected by using the Compile Notebook button, which runs where a script lives
## as opposed to the true home directory for this project
str(gCoef <- readRDS("../results/gCoef.rds"))
hDat <-
  droplevels(subset(gCoef,
                    continent %in% c("Asia", "Americas")))
str(hDat)
dotplot(slope ~ continent, hDat)
t.test(slope ~ continent, hDat)

## old school of doing something similar is this:
## at start place a "sink", like so:
## sink("slopeComparisonAsiaVsAmericas_fromSink.txt")
## <insert all the comands above here>
## sink()
## the file left behind is a (very) poor's man dynamic report