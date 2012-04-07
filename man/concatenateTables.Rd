\name{concatenateTables}
\alias{concatenateTables}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
concatenateTables(tables, timeOffsets = round(cumsum(sapply(tables, function(t) {
    return(aveVectorDelta(t$time) + tail(t$time, 1))
}))))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tables}{
%%     ~~Describe \code{tables} here~~
}
  \item{timeOffsets}{
%%     ~~Describe \code{timeOffsets} here~~
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
function (tables, timeOffsets = round(cumsum(sapply(tables, function(t) {
    return(aveVectorDelta(t$time) + tail(t$time, 1))
  })))) 
{
    l <- length(tables)
    fixedTables <- mapply(function(table, offset) {
        table$time <- table$time + offset
        return(table)
    }, tables[2:l], timeOffsets[1:(l - 1)], SIMPLIFY = FALSE)
    return(do.call(rbind, c(tables[1], fixedTables)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
