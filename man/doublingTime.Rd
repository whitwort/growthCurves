\name{doublingTime}
\alias{doublingTime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
doublingTime(table, filters = c(default.lagFilter, default.plateauFilter), number.format = round)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{table}{
%%     ~~Describe \code{table} here~~
}
  \item{filters}{
%%     ~~Describe \code{filters} here~~
}
  \item{number.format}{
%%     ~~Describe \code{number.format} here~~
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
function (table, filters = c(default.lagFilter, default.plateauFilter), 
    number.format = round) 
{
    calculateDT <- function(col) {
        filter <- c(TRUE)
        for (f in filters) {
            filter <- filter & f(col)
        }
        od <- col[filter]
        time <- table$time[filter]
        if (length(od) < 1) {
            return(NA)
        }
        fit <- lm(log(od) ~ time)
        return(number.format(log(2)/coef(fit)["time"]/60))
    }
    return(apply(table[2:length(table)], MARGIN = 2, FUN = calculateDT))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
