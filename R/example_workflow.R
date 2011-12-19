#-----------------------------------------------------------------------
# Example data workflow preprocessing in R
#
# Run this in the R/ directory of the pipeline dir
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# Libraries, Data and sources
#-----------------------------------------------------------------------

library(ggplot2)
library(KCsmart)

data(mmMirrorLocs)

source('functions_preprocessing.R')

#-----------------------------------------------------------------------
# Precomputation of variables
#-----------------------------------------------------------------------

# Cumulative chromosome end coordinates
chromEnds <- cumsum(unlist(lapply(mmMirrorLocs, max)))
chromEnds <- c(0, chromEnds)

#-----------------------------------------------------------------------
# Preprocessing main workflow
#-----------------------------------------------------------------------

# Read data
dataFile <- '../exampleData/24-06-09-3.export.txt'

# Import and filter/annotate the readFrame.
readFrame <- importReadFrame(dataFile)

# Do the clustering
readFrame <- readFrameInsertionClustering(readFrame)

# Remove contaminations
readFrame <- readFrameResolveContaminations(readFrame)

# Aggregate results into a custInfo frame
clustInfo <- readFrameFillClustInfo(readFrame)

# Write final result to file
write.table(clustInfo, file = "../exampleData/insInfo.csv", quote = FALSE, sep = ",", row.names = FALSE)

# Produce tumor information (reads per tumor)

tumorInfo <- data.frame(tumorIndex=vector(mode='character', length=unique(readFrame$index_numeric)), tumorIndexNum=unique(readFrame$index_numeric),
	totalReads=vector(mode='character', length=unique(readFrame$index_numeric)), stringsAsFactors=F)

for (i in 1:nrow(tumorInfo)) {
	tumorInfo$tumorIndex[i] <- unique(readFrame$index[readFrame$index_numeric == tumorInfo$tumorIndexNum[i]])
	tumorInfo$totalReads[i] <- sum(readFrame$depth[readFrame$index_numeric == tumorInfo$tumorIndexNum[i]])
}

write.table(tumorInfo, file = "../exampleData/tumInfo.csv", quote = FALSE, sep = ",", row.names = FALSE)
