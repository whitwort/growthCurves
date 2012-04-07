\name{magellanTable}
\alias{magellanTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
magellanTable(filePath, wellLabels = c(t(default.plate.96)), table.postprocess = validateTable)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filePath}{
%%     ~~Describe \code{filePath} here~~
}
  \item{wellLabels}{
%%     ~~Describe \code{wellLabels} here~~
}
  \item{table.postprocess}{
%%     ~~Describe \code{table.postprocess} here~~
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
function (filePath, wellLabels = c(t(default.plate.96)), table.postprocess = validateTable) 
{
    table <- read.table(filePath, sep = "\t", header = FALSE, 
        col.names = c("time", wellLabels, "RM"))
    ctime <- as.character(table$time)
    ntime <- as.numeric(substr(ctime, start = 1, stop = nchar(ctime) - 
        1))
    table$time <- ntime
    table$RM <- NULL
    return(table.postprocess(table))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
