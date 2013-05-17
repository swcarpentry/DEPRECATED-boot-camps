## this is an ugly script -- not meant to be "source'd" -- where I am merely
## collecting code we ran to learn about data aggregation
gDat <- read.delim("gapminderDataFiveYear.txt")
str(gDat)

## don't store copies of little bits of your data.frame like this unless you
## have good reason
(snippet <- subset(gDat, country == "Cambodia"))

## do you want to make plots for each value of a factor? then use multi-panel
## conditioning in lattice or facetting in ggplot2
library(lattice)
xyplot(lifeExp ~ year | country, gDat,
       subset = continent == "Oceania")

## let's begin with the data aggregation functions built-in to R

## let's begin with a function that operates on a matrix (or higher dimensional
## arrays, actually)
(jCountries <- sort(c("Canada", "Cambodia", "Rwanda")))
tinyDat <- subset(gDat, country %in% jCountries)
## JB inserted this later ... safety first! making sure jCountries is in same
## order as the countries appear in tinyDat (and, therefore, gDat)
## important when we apply column names below
(jCountries <- as.character(sort(unique(tinyDat$country))))
tinyDat <- matrix(tinyDat$lifeExp, ncol = length(jCountries))
colnames(tinyDat) <- jCountries
rownames(tinyDat) <- sort(unique(gDat$year))
apply(tinyDat, 1, mean)
apply(tinyDat, 2, median)
## the sensible names appearing here are the payoff for setting up colnames and
## rownames
rowMeans(tinyDat) # FAST use for big datasets
which.min(tinyDat[1, ])
jCountries[apply(tinyDat, 1, which.min)]

## moving on to sapply and lapply, which operate on lists
## recall that data.frames are a special case of a list
sapply(gDat, is.numeric)
gDatNum <- subset(gDat, select = sapply(gDat, is.numeric))
str(gDatNum)
head(gDatNum)
sapply(gDatNum, median)
lapply(gDatNum, median)
## notice how different the output looks
## why is that?
str(sapply(gDatNum, median))
str(lapply(gDatNum, median))
## lapply always returns a list
## sapply attempt to tidy up

## what if we want a summary that's 2-dimensional?
sapply(gDatNum, range)
lapply(gDatNum, range)

## moving to the case where we want to compute summaries on sub-data.frames,
## i.e. induced by the unique values of a factor (or combinations of multiple
## factors)
tapply(gDat$lifeExp, gDat$continent, max)
with(gDat, tapply(lifeExp, continent, max))
with(gDat,
     tapply(country, continent, function(x) {
       length(unique(x))
     }))

## the un-predictability of what apply, sapply, lapply, tapply can return is
## disappointing

## for data aggregation BLISS, use the add-on package plyr

## JB finds hard to learn from documentation of individual functions
## read this paper to get the big picture:
## http://www.jstatsoft.org/v40/i01/paper

install.packages(pkgs = "plyr")
library(plyr)

ddply(gDat, .(continent), summarise, median = median(lifeExp))

str(ddply(gDat, .(continent), summarise, median = median(lifeExp)))
str(tapply(gDat$lifeExp, gDat$continent, median))


ddply(gDat, .(continent), summarise,
      min = min(lifeExp), max = max(lifeExp))

## let's run linear regression of lifeExp on year for individual countries and
## save the estimated intercept and slope

## walk before you run ....
lm(lifeExp ~ year, gDat, subset = country == "Zimbabwe")
xyplot(lifeExp ~ year, gDat, subset = country == "Zimbabwe")
yearMin <- min(gDat$year)
lm(lifeExp ~ I(year - yearMin), gDat,
   subset = country == "Zimbabwe")
coef(lm(lifeExp ~ I(year - yearMin), gDat,
        subset = country == "Zimbabwe"))

## package your working protoype code in a function
jFun <- function(z) {
  yearMin <- min(z$year)
  jCoef <- coef(lm(lifeExp ~ I(year - yearMin), z))
  names(jCoef) <- c("intercept", "slope")
  return(jCoef)
}

## test your function!
jFun(subset(gDat, country == "Canada"))

## scale up:
## let ddply to handle all the booking keeping, i.e. managing the loop over
## countries
gCoef <- ddply(gDat, .(country), jFun)
str(gCoef)
tail(gCoef)

## I wish that I also had the continent info
## two ways to get that:

## easiest:
## sort of a trick: add continent to the ddply call
gCoef <- ddply(gDat, .(country, continent), jFun)
str(gCoef)
tail(gCoef)

## use match()
gCoef <- ddply(gDat, .(country), jFun)
str(gCoef)
tail(gCoef)
gCoef$continent <- gDat$continent[match(gCoef$country, gDat$country)]
str(gCoef)
tail(gCoef)
