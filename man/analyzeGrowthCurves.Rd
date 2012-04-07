\name{analyzeGrowthCurves}
\alias{analyzeGrowthCurves}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
analyzeGrowthCurves(tablePath, savePath = NULL, tableParser = magellanTable, annotationPath = NULL, annotationParser = plateSetup, plateLabels = default.plate.96, filters = c(default.lagFilter, default.plateauFilter), zipPath = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tablePath}{
%%     ~~Describe \code{tablePath} here~~
}
  \item{savePath}{
%%     ~~Describe \code{savePath} here~~
}
  \item{tableParser}{
%%     ~~Describe \code{tableParser} here~~
}
  \item{annotationPath}{
%%     ~~Describe \code{annotationPath} here~~
}
  \item{annotationParser}{
%%     ~~Describe \code{annotationParser} here~~
}
  \item{plateLabels}{
%%     ~~Describe \code{plateLabels} here~~
}
  \item{filters}{
%%     ~~Describe \code{filters} here~~
}
  \item{zipPath}{
%%     ~~Describe \code{zipPath} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (tablePath, savePath = NULL, tableParser = magellanTable, 
    annotationPath = NULL, annotationParser = plateSetup, plateLabels = default.plate.96, 
    filters = c(default.lagFilter, default.plateauFilter), zipPath = NULL, 
    ...) 
{
    optArgs <- list(...)
    optArgs$wellLabels = c(t(plateLabels))
    wrapCall <- function(f) {
        return(cleanWrapper(f, optArgs))
    }
    filters <- lapply(filters, wrapCall)
    optArgs$filters <- filters
    tables <- lapply(tablePath, wrapCall(tableParser))
    if (length(tables) > 1) {
        table <- concatenateTables(tables)
    }
    else {
        table <- tables[[1]]
    }
    if (!is.null(annotationPath)) {
        resultTable <- wrapCall(annotationParser)(filePath = annotationPath)
    }
    else {
        resultTable <- data.frame(row.names = optArgs$wellLabels, 
            check.names = TRUE)
    }
    resultTable$doublingTime <- wrapCall(doublingTime)(table = table)
    optArgs$annotations <- resultTable
    if (!is.null(savePath)) {
        wrapCall(makeODPlots)(table = table, savePath = savePath)
    }
    if (!is.null(savePath)) {
        write.table(resultTable, paste(savePath, "results.txt", 
            sep = "/"), sep = "\t", col.names = NA)
    }
    if (!is.null(zipPath) && !is.null(savePath)) {
        zip(zipPath, files = paste(savePath, dir(savePath), sep = "/"))
    }
    return(resultTable)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
