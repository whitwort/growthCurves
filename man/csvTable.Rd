\name{csvTable}
\alias{csvTable}
\title{Loads data from a csv (comma-separated) text file.}
\usage{
csvTable(filePath)
}
\arguments{
  \item{filePath}{path on the local file system to be
  passed to the parser}
}
\value{
a data.frame object holding the parsed data.
}
\description{
This function uses the builtin \code{\link{read.csv}} to
parse a comma-separated file containg OD density readings
over the experimental time course.  The left- most column
should contain time points (with earliest-to-latest
timepoints givendown the rows of the table).  The numbers
in this column should be unquoted. The remaining columns
should give, in order, the time course data for each well
on the plate (for example, starting with well A1 in column
to and ending with well H12 in column 97).
}

