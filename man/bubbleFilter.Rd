\name{bubbleFilter}
\alias{bubbleFilter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bubbleFilter(od, bubble.tolerance = 3, bubble.neighborhoodSize = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{od}{
%%     ~~Describe \code{od} here~~
}
  \item{bubble.tolerance}{
%%     ~~Describe \code{bubble.tolerance} here~~
}
  \item{bubble.neighborhoodSize}{
%%     ~~Describe \code{bubble.neighborhoodSize} here~~
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
function (od, bubble.tolerance = 3, bubble.neighborhoodSize = 3) 
{
    localDelta <- vectorDeltas(od)
    neighborhoodDelta <- vector("numeric", length(localDelta))
    for (i in 1:length(localDelta)) {
        neighborhoodDelta[i] <- abs(localDelta[i] - median(localDelta[max(1, 
            i - bubble.neighborhoodSize):min(length(localDelta), 
            i + bubble.neighborhoodSize)]))
    }
    cutoff <- mean(neighborhoodDelta) * (1 + bubble.tolerance)
    goodPoints <- (neighborhoodDelta < cutoff)
    return(goodPoints)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
