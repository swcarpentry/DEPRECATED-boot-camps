## this was our first complete small "analysis" (of simulated data)
## we learned how to take code we typed into the Console,
## select it in the History, and send it to Source, i.e. the
## script editor pane
## there we tweaked it a bit, saved to this file, and practiced
## using Run, Source, etc via mouse and keyboard shortcuts
n <- 100
a <- 2
b <- 3
sigSq <- 0.5
x <- runif(n)
y <- a + b * x + rnorm(n, sd = sqrt(sigSq))
(avgX <- mean(x))
write(avgX, "avgX.txt")
plot(x, y)
abline(a, b, col = "purple", lwd = 3)
dev.print(pdf, "niftyPlot.pdf")
