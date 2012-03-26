#growthCurves

growthCurves is a small [R script](http://www.r-project.org/) designed to aid analysis of liquid culture microbial growth curves obtained in 96- or 384-well format.

Deeper documentation will be comming soon, but what follows should be enough to get started if you are familiar with R.


##Loading the script
Although growthCurves is not yet setup as a full R package, it is currently implemented in a single R script so it is easy to load into your current R session using the 'source()' function.  If you'd like to stay uptodate with the latest development version you can include the something like the following in your code:

```r
require(RCurl)
url <- "https://github.com/whitwort/growthCurves/blob/master/growthCurves.R"
eval(parse(text = getURL(url, ssl.verifypeer = FALSE)), envir=.GlobalEnv)
```

Unfortunately, using source(url) directly won't work because 'source' doesn't support HTTPS used by github).  You'll need to have the CRAN package 'RCurl' installed for this to work.


##Running an analysis
The wrapper function 'analyzeGrowthCurves' provides a simple interface to running a full analysis which produces growth curve graphs, incorporates well annotations, and calculates a doubling time for each well.

With the [example files](https://github.com/whitwort/growthCurves/tree/master/examples) included in this repository in the current working directory, this function could be called with:

```r
analyzeGrowthCurves( tablePath="sample.asc", annotationPath="plateSetup.txt", savePath="" )
```

This will run an analysis producing an number of files in the current working directory.  'results.txt' will hold a tab-delimited table incorporating well annotations and the calculated doubling times.  Images files will also be produced showing each growth curve individually or as a grid in 'composite.png'.


## License

Copyright © 2012 Gregg Whitworth.  Licensed under the [MIT License](http://mit-license.org).