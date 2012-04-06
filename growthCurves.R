#Wrapper functions
analyzeGrowthCurves <- function(tablePath, 
                                savePath          = NULL,
                                tableParser       = magellanTable, 
                                annotationPath    = NULL, 
                                annotationParser  = plateSetup,
                                plateLabels       = default.plate.96,
                                filters           = c(default.lagFilter, default.plateauFilter),
                                zipPath           = NULL,
                                ...
                                )  {
  
  #Listify additional optional arguments
  optArgs <- list(...)
  
  #Calculate common arguments
  optArgs$wellLabels = c(t(plateLabels))

  #We'll simplify the wrapping interface with a little closure
  wrapCall <- function(f) { return( cleanWrapper(f,optArgs) ) }
  
  #Automatically wrap the filters, and save to optArgs
  filters <- lapply(filters, wrapCall)
  optArgs$filters <- filters
  
  #Load the OD table(s) with the given parser
  tables <- lapply(tablePath, wrapCall(tableParser))
  
  #If there is more than one table, concatenate
  if (length(tables) > 1) {
    table <- concatenateTables(tables)
  } else {
    table <- tables[[1]]
  }
  
  #Load the annotations if we're given a file path
  if (!is.null(annotationPath)) {
    resultTable <- wrapCall(annotationParser)(filePath = annotationPath)
  } else {
    resultTable <- data.frame(row.names=optArgs$wellLabels, check.names = TRUE)
  }
  
  #Run analysis functions
  resultTable$doublingTime <- wrapCall(doublingTime)(table = table)
  
  #Save analysis and original annotations to optArgs
  optArgs$annotations <- resultTable
  
  #Run visualization functions, if we were given a save path
  if (!is.null(savePath)) {
    wrapCall(makeODPlots)(table = table, savePath = savePath)
  }
  
  #Save the analysis and annotations out to a result table and return it
  if (!is.null(savePath)) {
    write.table(resultTable, paste(savePath, "results.txt", sep="/"), sep="\t", col.names=NA)
  }
  
  #If we were given a zip path
  if (!is.null(zipPath) && !is.null(savePath)) {
    zip( zipPath, files = dir(savePath) )
  }
  
  return(resultTable)
  
}

#Default labels for 96- and 384-well plates (character matrix with the same shape as the plate)
default.plate.96 = matrix(
  data  = paste( rep(LETTERS[1:8], each = 12), 1:12, sep = ""), 
  nrow  = 8, 
  ncol  = 12, 
  byrow = TRUE
  )

default.plate.384 = matrix(
  data  = paste( rep(LETTERS[1:16], each = 24), 1:24, sep = ""), 
  nrow  = 16, 
  ncol  = 24, 
  byrow = TRUE
  )

#OD table parsers
csvTable <- function(filePath) { return( read.csv(filePath) ) }
tabTable <- function(filePath) { return( read.delim(filePath) ) }
xlsTable <- function(filePath) {
  #This function requires that both gdata and 'perl' be installed
  library(gdata)
  return( read.xls(filePath) )
}

#TODO Magellan export instructions
magellanTable <- function(filePath, 
                          wellLabels        = c(t(default.plate.96)),
                          table.postprocess = validateTable
                          ) {
  
  #Read the table, using labels as the column names
  table <- read.table(filePath, 
                      sep       = "\t", 
                      header    = FALSE, 
                      col.names = c("time", wellLabels, "RM"))
  
  #Convert the first column to a numeric vector by stripping the "s" from each entry
  ctime <- as.character(table$time)
  ntime <- as.numeric(substr(ctime, start = 1, stop = nchar(ctime)-1 ))
  
  table$time <- ntime
  
  #Cleanup the spurious last column
  table$RM <- NULL
  
  return(table.postprocess(table))
}

validateTable <- function(table) {
  
  #Check each row to see if it contains NA's
  badRows <- apply(table, 1, function(r) { NA %in% r })
  
  #If there are bad rows
  if (TRUE %in% badRows) {
    print("The following timepoints will be purged because they contain incomplete data:")
    print(paste(table$time[badRows], sep = ","))
    
    #Dump only good rows back into table
    table <- table[!badRows,]
  }
  
  return(table)
  
}

concatenateTables <- function(tables, 
                              timeOffsets = round( cumsum( sapply( tables, 
                                            function(t) { 
                                              return(
                                                aveVectorDelta(t$time) + tail(t$time,1)
                                                )}
                                            )))
                              ) {
  
  #Apply each offset to the times on each table
  l <- length(tables)
  fixedTables <- mapply(
    function(table, offset) {
      table$time <- table$time + offset
      return(table)
    },
    
    #Apply first offset to the second table, and so on...
    tables[2:l], 
    timeOffsets[1:(l-1)], 
    SIMPLIFY = FALSE
    
  )
  
  #Concatenate the data frames using rbind (row bind), executed as a do.call
  return( do.call(rbind, c(tables[1], fixedTables)) )
  
}

#Plate annotation parsers
plateSetup <- function(filePath, plateLabels = default.plate.96) {
  
  #Parsers
  parseBlock <- function(block) {
    
    #Parse each line in the block, return a named vector
    return(unlist(lapply(block, parseLine)))
    
  }
  
  parseLine <- function(line) {
    
    #Extract the well definition text and the category label from this line
    m <- strsplit(line, "(\\s*):(\\s*)")[[1]]
    wells <- parseWells(m[1])
    labels <- rep(m[2], times = length(wells))
    
    #Assign labels to the wells vector as names
    names(labels) <- wells
    
    #Return a labels vector named by associated well
    return(labels)
    
  }
  
  parseWells <- function(wells) {
    
    #If we're passed a single well label return it
    if (wells %in% plateLabels) { return(wells) }
    
    #If stripping all whitespace gives us a label, return it
    stripped <- gsub("\\s", "", wells)
    if (stripped %in% plateLabels) { return( stripped ) }
    
    #If we're passed a block range ("->")
    rangeSplit <- strsplit(wells, "(\\s*)->(\\s*)")[[1]]
    if (length(rangeSplit) == 2) {
      
      #Get the start and ending well labels
      rangeStart <- parseWells(rangeSplit[1])
      rangeEnd <- parseWells(rangeSplit[2])
      
      #Convert these labels to indexes in the matrix
      startIndex <- which(plateLabels == rangeStart, arr.ind=TRUE)
      endIndex <- which(plateLabels == rangeEnd, arr.ind=TRUE)
      
      #Return the block of labels between startIndex and endIndex
      return(  
        as.vector(
          plateLabels[startIndex[1]:endIndex[1],startIndex[2]:endIndex[2]]
          )
        )
    }
    
    #TODO handle or's
    print(paste("Could not parse well:", wells))
    
  }

  
  #Load the file by line, ignoring lines with only whitespace
  fileLines <- scan(filePath, character(0), sep = "\n", strip.white=TRUE)
  
  #Calculate a vector containing the indexes of each block heading in lines 
  #where blocks headings are defined by lines without a ":"
  blockHeadings <- (1:length(fileLines))[!grepl(":", fileLines, fixed=TRUE)]
  
  #For each block heading...
  annotations <- data.frame(row.names = c(t(plateLabels)), check.rows = TRUE)
  for (i in 1:length(blockHeadings)) {
    
    #Get this and the next heading index
    thisHeadingIndex <- blockHeadings[i]
    nextHeadingIndex <- min(blockHeadings[i+1], length(fileLines)+1, na.rm = TRUE)
    blockName <- fileLines[thisHeadingIndex]
    
    #Parse this block
    wellAnnotations <- parseBlock(fileLines[(thisHeadingIndex+1):(nextHeadingIndex-1)])
    
    #Check to make sure each well is defined only once
    if ( length(names(wellAnnotations)) != length(unique(names(wellAnnotations))) ) {
      print(paste( 
        "Warning:  in block '",  blockName, "' there are wells which are given more than one label."
        ))  
    }
    
    #Add the annotation set to the list of annotations, reordering wellLabels
    annotations[[blockName]] <- wellAnnotations[c(t(plateLabels))]
    
  }
  
  return(annotations)
  
}

#Filter Functions, return logical vectors where FALSE indicates bad data
default.lagFilter <- function(od, lag.window = 3) {
  return(od > (median(od[1:lag.window]) * 2))
}

default.plateauFilter <- function(od, plateau.cutoff = 0.85) {
  return(od < plateau.cutoff)
}

bubbleFilter <- function(od, bubble.tolerance = 3, bubble.neighborhoodSize = 3) {
  
  #Calculate pair-wise point-to-point value changes (deltas)
  localDelta <- vectorDeltas(od)
  
  #Compare each pair-wise delta to the average delta in its neighborhood
  neighborhoodDelta <- vector("numeric", length(localDelta))
  for (i in 1:length(localDelta)) {
    neighborhoodDelta[i] <- abs( 
      localDelta[i] - median(
        localDelta[max(1,i-bubble.neighborhoodSize):min(length(localDelta),i+bubble.neighborhoodSize)])
      )
  }
  
  #Caclulate a neighborhood cutoff value based on the given tolerance parameter
  cutoff <- mean(neighborhoodDelta) * (1 + bubble.tolerance)
  
  #Keep only points with neighborhood deltas below the cutoff
  goodPoints <- (neighborhoodDelta < cutoff)
  
  #Return the final set of points passing the filter
  return( goodPoints )
  
}

#Analysis functions
doublingTime <- function(table, 
                         filters        = c(default.lagFilter, default.plateauFilter), 
                         number.format  = round
                         ) {
  
  calculateDT <- function(col) {
     #Apply each filter to the column of OD data
    filter <- c(TRUE)
    for (f in filters) { filter <- filter & f(col) }
    
    #Use ODs and times that pass the filters
    od <- col[filter]
    time <- table$time[filter]

    #Check to ensure we haven't filtered all points
    if (length(od)<1) { return(NA) }
    
    #Fit a linear model to the log transformed ODs versus time
    fit <- lm(log(od) ~ time)
    
    #Calculate the doubling time
    return( number.format(log(2) / coef(fit)["time"] / 60 ) )
  }
  
  return( apply(table[2:length(table)], MARGIN = 2, FUN = calculateDT) )
  
}

#Visualizations
makeODPlots <- function(table, savePath, 
                        annotations = data.frame(),
                        plateLabels = default.plate.96, 
                        filters     = c(default.lagFilter, default.plateauFilter),
                        individual  = TRUE,
                        composite   = c(300,350)
                        ) {
  
  #Calculate a linear series of labels
  wellLabels <- c(t(plateLabels))
  
  #Calculate a universal set of x- and y- axis scales
  xlim <- c( min(table$time), max(table$time) ) / 60
  ylim <- c( min(table[2:length(table)]), max(table[2:length(table)]) )
  
  #od Plotting function (used to make both individual and aggregated graphs)
  odPlot <- function(well, 
                     main     = NULL, 
                     xlab     = NULL, 
                     ylab     = NULL, 
                     ul.label = NULL, 
                     lr.label = NULL
                     ) {
    
    #Apply each filter to the column of OD data
    filter <- c(TRUE)
    for (f in filters) { filter <- filter & f(table[[well]]) }
    
    #Plot the good points in blue
    convertedTimes <- table$time[filter] / 60
    plot(convertedTimes, table[[well]][filter], 
         main = main,
         xlab = xlab, 
         ylab = ylab, 
         xlim = xlim, 
         ylim = ylim,
         col  = 4,
         pch  = 20,
         )
    
    #Plot the filtered points in black
    points(table$time[!filter] / 60, table[[well]][!filter], pch = 20, col = "grey45")
    
    #If a doubling time is saved on annotations, plot the doubling time function
    if (!is.null(annotations$doublingTime)) {
      dt <- annotations[well,'doublingTime']
      if (!is.na(dt)) {
        
        startOD <- min(table[[well]][filter])
        logStartTime <- convertedTimes[1]
        
        model <- function(t) { startOD * ( 2^(t/dt) ) }
        lines(table$time/60, model((table$time/60)-logStartTime), col=4)
        
      }
    }
    
    #If we were given an upper left annotation label
    if (!is.null(ul.label)) {  
      #Put the annotation text in the upper lefthand corner of the graph
      text(x = xlim[1], y = ylim[2], labels = ul.label, adj = c(0,1))
    }
    
    #If we were given an upper left annotation label
    if (!is.null(lr.label)) {
      #Put the annotation text in the lower righthand corner of the graph
      text(x = xlim[2], y = ylim[1], labels = ul.label, adj = c(1,0))
    }
    
  }
  
  #Make a plot of each individual OD series versus time
  if (individual) {
    for (well in wellLabels) {
      
      #If we were passed annotations, write them to the graph as the upper left label
      if (length(annotations)>0) {
        ul.label <- paste(
          names(annotations), ": ", annotations[well,], sep = "", collapse = "\n" 
          )
      } else { 
        ul.label = NULL
      }
      
      #Start a new JPEG
      jpeg(paste(savePath, paste(well, ".jpg"), sep="/"))
      
      #Make the individual plot
      odPlot(well, main = well, xlab = "Time (minutes)", ylab = "OD", ul.label = ul.label)
      
      #Save the image
      dev.off()
      
    }
  }
  
  #Make a composite plot of all series
  if (composite[1]) {
    
    #Calculate a size for the composite image and open a file for writing
    imageDim <- composite * dim(plateLabels) #c(row, col)
    png(filename=paste(savePath, "composite.png", sep = "/"), width = imageDim[2], height = imageDim[1])
    
    #Save the starting par values, so that we can restore later
    original.par = par(
      mfrow   = dim(plateLabels),
      mar     = c(0.3,0.3,0.3,0.3),
      oma     = c(4,4,4,4), 
      yaxt    = "n",
      xaxt    = "n",
      cex     = 1
      )
    
    #For each plot
    for (well in wellLabels) {
      
      #If we were passed annotations, write them to the graph as the upper left label
      if (length(annotations)>0) {
        ul.label <- paste(c(well, annotations[well,]), sep = "", collapse = "\n")
      } else { 
        ul.label = NULL
      }
      
      #Make the plot
      odPlot(well, xlab="", ylab="", ul.label=ul.label)
      
    }
    
    #Save the image
    dev.off()
    
    #Restore the original plot parameters
    par(original.par)
    
  }
    
}

#Utility functions
vectorDeltas <- function(v) { return( abs( v[-1] - v[-length(v)] ) ) }
aveVectorDelta <- function(v, ave = mean) { return( ave(vectorDeltas(v)) ) }
cleanCall <- function(f, ..., optArgs=list()) {
  #TODO is there anything in the standard library like this?
  
  #If there are optional arguments, pass through those that match to f's arguments
  if (!is.null(names(optArgs))) {
    goodArgs <- !is.na(pmatch( names(optArgs), names(formals(f)) ))
    callList <- c(list(...), optArgs[goodArgs])
  } else {
    callList <- list(...)
  }
  
  #Call with the good arguments, and pass-through the result value
  return( do.call(f, callList) )
  
}
cleanWrapper <- function(f, optArgs) {
  force(f)
  return(
    function(...) {
      return( cleanCall(f, ..., optArgs = optArgs) )
    }
    )
}
