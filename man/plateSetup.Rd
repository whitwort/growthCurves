\name{plateSetup}
\alias{plateSetup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plateSetup(filePath, plateLabels = default.plate.96)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filePath}{
%%     ~~Describe \code{filePath} here~~
}
  \item{plateLabels}{
%%     ~~Describe \code{plateLabels} here~~
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
function (filePath, plateLabels = default.plate.96) 
{
    parseBlock <- function(block) {
        return(unlist(lapply(block, parseLine)))
    }
    parseLine <- function(line) {
        m <- strsplit(line, "(\\s*):(\\s*)")[[1]]
        wells <- parseWells(m[1])
        labels <- rep(m[2], times = length(wells))
        names(labels) <- wells
        return(labels)
    }
    parseWells <- function(wells) {
        if (wells \%in\% plateLabels) {
            return(wells)
        }
        stripped <- gsub("\\s", "", wells)
        if (stripped \%in\% plateLabels) {
            return(stripped)
        }
        rangeSplit <- strsplit(wells, "(\\s*)->(\\s*)")[[1]]
        if (length(rangeSplit) == 2) {
            rangeStart <- parseWells(rangeSplit[1])
            rangeEnd <- parseWells(rangeSplit[2])
            startIndex <- which(plateLabels == rangeStart, arr.ind = TRUE)
            endIndex <- which(plateLabels == rangeEnd, arr.ind = TRUE)
            return(as.vector(plateLabels[startIndex[1]:endIndex[1], 
                startIndex[2]:endIndex[2]]))
        }
        print(paste("Could not parse well:", wells))
    }
    fileLines <- scan(filePath, character(0), sep = "\n", strip.white = TRUE)
    blockHeadings <- (1:length(fileLines))[!grepl(":", fileLines, 
        fixed = TRUE)]
    annotations <- data.frame(row.names = c(t(plateLabels)), 
        check.rows = TRUE)
    for (i in 1:length(blockHeadings)) {
        thisHeadingIndex <- blockHeadings[i]
        nextHeadingIndex <- min(blockHeadings[i + 1], length(fileLines) + 
            1, na.rm = TRUE)
        blockName <- fileLines[thisHeadingIndex]
        wellAnnotations <- parseBlock(fileLines[(thisHeadingIndex + 
            1):(nextHeadingIndex - 1)])
        if (length(names(wellAnnotations)) != length(unique(names(wellAnnotations)))) {
            print(paste("Warning:  in block '", blockName, "' there are wells which are given more than one label."))
        }
        annotations[[blockName]] <- wellAnnotations[c(t(plateLabels))]
    }
    return(annotations)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
