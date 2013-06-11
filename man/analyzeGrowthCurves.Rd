\name{analyzeGrowthCurves}
\alias{analyzeGrowthCurves}
\title{Kitchen-sink wrapper for optical growth curve analysis.}
\usage{
  analyzeGrowthCurves(tablePath, savePath = NULL,
    tableParser = magellanTable, annotationPath = NULL,
    annotationParser = plateSetup,
    plateLabels = default.plate.96,
    filters = c(default.lagFilter, default.plateauFilter),
    zipPath = NULL, ...)
}
\arguments{
  \item{tablePath}{path on the local file system to be
  passed to the tableParser}

  \item{savePath}{optional path to a directory where
  results should be saved; defaults to NULL which means
  that no tables or images will be saved.}

  \item{tableParser}{a function to use to parse tables; see
  \code{\link{csvTable}} for an example.}

  \item{annotationPath}{path on the local file system to
  pass to the annotationParser}

  \item{annotationParser}{a function to use to parse
  annotation files; see \code{\link{plateSetup}} for an
  example.}

  \item{plateLabels}{a matrix of well labels that matches
  the shape of the plate data being analyzed.  Two default
  layouts are provided: \code{default.plate.96} and
  \code{\link{default.plate.384}}}

  \item{filters}{a list of unnamed functions which will be
  used to filter OD data; see
  \code{\link{default.lagFilter}} and
  \code{\link{default.plateauFilter}} for examples.}

  \item{zipPath}{path on the local file system where a new
  zipfile archive should be written containing all files
  produced in the analysis.  Defaults to NULL which
  indicates no zip file should be created.}

  \item{...}{all additional arguments are passed along to
  the parser, filter, analysis and visualization functions
  that recieve them.}
}
\value{
  a data.frame object holding the results of the analysis
  and annotations
}
\description{
  This function loads one or more tables, a set of well
  label annotations (if given), and then calculates
  doubling times (\code{\link{doublingTime}}) and produces
  individual and composite OD plots
  (\code{\link{makeODPlots}}). The default analysis
  pipeline can be tuned by passing any arguments accepted
  by table parsers, annotation parsers, data filters,
  analysis functions and visualization functions.
}

