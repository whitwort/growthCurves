\name{makeODPlots}
\alias{makeODPlots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makeODPlots(table, savePath, annotations = data.frame(), plateLabels = default.plate.96, filters = c(default.lagFilter, default.plateauFilter), individual = TRUE, composite = c(300, 350))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{table}{
%%     ~~Describe \code{table} here~~
}
  \item{savePath}{
%%     ~~Describe \code{savePath} here~~
}
  \item{annotations}{
%%     ~~Describe \code{annotations} here~~
}
  \item{plateLabels}{
%%     ~~Describe \code{plateLabels} here~~
}
  \item{filters}{
%%     ~~Describe \code{filters} here~~
}
  \item{individual}{
%%     ~~Describe \code{individual} here~~
}
  \item{composite}{
%%     ~~Describe \code{composite} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (table, savePath, annotations = data.frame(), plateLabels = default.plate.96, 
    filters = c(default.lagFilter, default.plateauFilter), individual = TRUE, 
    composite = c(300, 350)) 
{
    wellLabels <- c(t(plateLabels))
    xlim <- c(min(table$time), max(table$time))/60
    ylim <- c(min(table[2:length(table)]), max(table[2:length(table)]))
    odPlot <- function(well, main = NULL, xlab = NULL, ylab = NULL, 
        ul.label = NULL, lr.label = NULL) {
        filter <- c(TRUE)
        for (f in filters) {
            filter <- filter & f(table[[well]])
        }
        convertedTimes <- table$time[filter]/60
        plot(convertedTimes, table[[well]][filter], main = main, 
            xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, 
            col = 4, pch = 20, )
        points(table$time[!filter]/60, table[[well]][!filter], 
            pch = 20, col = "grey45")
        if (!is.null(annotations$doublingTime)) {
            dt <- annotations[well, "doublingTime"]
            if (!is.na(dt)) {
                startOD <- min(table[[well]][filter])
                logStartTime <- convertedTimes[1]
                model <- function(t) {
                  startOD * (2^(t/dt))
                }
                lines(table$time/60, model((table$time/60) - 
                  logStartTime), col = 4)
            }
        }
        if (!is.null(ul.label)) {
            text(x = xlim[1], y = ylim[2], labels = ul.label, 
                adj = c(0, 1))
        }
        if (!is.null(lr.label)) {
            text(x = xlim[2], y = ylim[1], labels = ul.label, 
                adj = c(1, 0))
        }
    }
    if (individual) {
        for (well in wellLabels) {
            if (length(annotations) > 0) {
                ul.label <- paste(names(annotations), ": ", annotations[well, 
                  ], sep = "", collapse = "\n")
            }
            else {
                ul.label = NULL
            }
            jpeg(paste(savePath, paste(well, ".jpg", sep = ""), 
                sep = "/"))
            odPlot(well, main = well, xlab = "Time (minutes)", 
                ylab = "OD", ul.label = ul.label)
            dev.off()
        }
    }
    if (composite[1]) {
        imageDim <- composite * dim(plateLabels)
        png(filename = paste(savePath, "composite.png", sep = "/"), 
            width = imageDim[2], height = imageDim[1])
        original.par = par(mfrow = dim(plateLabels), mar = c(0.3, 
            0.3, 0.3, 0.3), oma = c(4, 4, 4, 4), yaxt = "n", 
            xaxt = "n", cex = 1)
        for (well in wellLabels) {
            if (length(annotations) > 0) {
                ul.label <- paste(c(well, annotations[well, ]), 
                  sep = "", collapse = "\n")
            }
            else {
                ul.label = NULL
            }
            odPlot(well, xlab = "", ylab = "", ul.label = ul.label)
        }
        dev.off()
        par(original.par)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
