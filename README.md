#growthCurves

A small [R script](http://www.r-project.org/) designed to aid analysis of liquid culture microbial growth curves obtained in 96- or 384-well format.

Deeper documentation will be comming soon, but what follows should be enough to get started if you are familiar with R.

##Loading the script

Although growthCurves.R is not yet available as a R package, it is implemented in a single R script which can be easily loaded with the 'source' function.  To load the latest development version of the script from this Github repository use:

```javascript
source(https://github.com/whitwort/growthCurves/blob/master/growthCurves.R)
```

##Running an analysis
The wrapper function 'analyzeGrowthCurves' provides a simple interface to running a full analysis which produces growth curve graphs, incorporates well annotations, and calculates a doubling time for each well.

With the [example files](/examples) included in this repository in the current working directory, this function could be called with:

```javascript
analyzeGrowthCurves( tablePath="sample.asc", annotationPath="plateSetup.txt", savePath="")
```

This will run an analysis producing an number of files in the current working directory.  'results.txt' will hold a tab-delimited table incorporating well annotations and the calculated doubling times.  Images files will also be produced showing each growth curve individually or as a grid in 'composite.png'.