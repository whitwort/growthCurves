#growthCurves

growthCurves is a small [R script](http://www.r-project.org/) designed to aid analysis of liquid culture microbial growth curves obtained in 96- or 384-well format.

Deeper documentation will be comming soon, but what follows should be enough to get started if you are familiar with R.

##Installing

To install the development version of this package run:

```r
devtools::install_github("whitwort/growthCurves")
```

##Running an analysis
The wrapper function 'analyzeGrowthCurves' provides a simple interface to running a full analysis which produces growth curve graphs, incorporates well annotations, and calculates a doubling time for each well.

With the [example files](https://github.com/whitwort/growthCurves/tree/master/examples) included in this repository in the current working directory, this function could be called with:

```r
library(growthCurves)

analyzeGrowthCurves( tablePath="sample.asc", annotationPath="plateSetup.txt", savePath="" )
```

This will run an analysis producing an number of files in the current working directory.  'results.txt' will hold a tab-delimited table incorporating well annotations and the calculated doubling times.  Image files will also be produced showing each growth curve individually or as a grid in 'composite.png'.


## License

Copyright © 2012 Gregg Whitworth and licensed under [LGPLv3](http://www.gnu.org/copyleft/lesser.html).