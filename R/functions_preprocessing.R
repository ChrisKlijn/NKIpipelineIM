#-----------------------------------------------------------------------
# Functions used in the R preprocessing pipeline for MMTV insertional
# mutagenesis screens using 454 sequencing.
#-----------------------------------------------------------------------

# TO DO: 
# 1. Start adding attributes to readFrame filtering to see how many
# reads were discarded caused by what
# 2. Make filtering by unique ligation point optional.

readFrameInsertionClustering <- function (readFrame) {

		# Separate into minus and plus strands
		readFrame_all_plus <- readFrame[readFrame$ori == '+',]
		readFrame_all_minus <- readFrame[readFrame$ori == '-',]

		# Initializa a cluster field for each insertion
		readFrame_all_plus$cluster <- 0
		readFrame_all_minus$cluster <- 0

		# Do the clustering for plus
		rf.dist <- as.dist(sapply(readFrame_all_plus$start_con, function (x) { abs(readFrame_all_plus$start_con - x) } ))
		rf.hc <- hclust(rf.dist, method='single')
    readFrame_all_plus$cluster <- cutree(rf.hc, h=10)

		# Get offset for the minus Cluster info
		maxClustPlus <- max(readFrame_all_plus$cluster)

		# Do the clustering for minus
		rf.dist <- as.dist(sapply(readFrame_all_minus$end_con, function (x) { abs(readFrame_all_minus$end_con - x) } ))
		rf.hc <- hclust(rf.dist, method='single')
    readFrame_all_minus$cluster <- cutree(rf.hc, h=10) + maxClustPlus

		readFrame_output <- rbind(readFrame_all_plus, readFrame_all_minus)
		return(readFrame_output)
}
	
readFrameFillClustInfo <- function (readFrame) {

	uniqClust <- unique(readFrame$cluster)
	clustNum <- length(uniqClust)
	
	# Initialize the clustinfo
	clustInfo <- data.frame(uniqClust=uniqClust, 
		uniqLig=vector(mode='numeric', length=clustNum), depth=vector(mode='numeric', length=clustNum), 
		stringsAsFactors=F, chr=vector(mode='numeric', length=clustNum), loc=vector(mode='numeric', length=clustNum),
		BC=vector(mode='numeric', length=clustNum), TumorID=vector(mode='character', length=clustNum),
		library=vector(mode='character', length=clustNum), strand=vector(mode='character', length=clustNum))
	
	for (i in 1:length(uniqClust)) {
		readFrame_clust <- readFrame[readFrame$cluster == uniqClust[i],]
    
    # If there is any read with an T7 flag, use the number of rows with the flag
    # to determine the Unique Ligation points, otherwise assign 1 to the UL
    if (any(readFrame_clust$flag == 1)) {
      clustInfo$uniqLig[i] <- nrow(readFrame_clust[readFrame_clust$flag == 1,])
    }
    else {
      clustInfo$uniqLig[i] <- 1
    }
		clustInfo$depth[i] <- sum(readFrame_clust$depth)
		clustInfo$chr[i] <- unique(readFrame_clust$numchr)
		clustInfo$BC[i] <- readFrame_clust$index_numeric[1]
		if (readFrame_clust$ori[1] == '+') {
			clustInfo$loc[i] <- readFrame_clust$start[1]
		}
		else {
			clustInfo$loc[i] <- readFrame_clust$end[1]
		}
		# If the data has been annotated with tumor information (eg.
		# unique tumor IDs and genotypes) add that to the clustInfo as
		# well
		if ('TumorID' %in% colnames(readFrame)) {
			clustInfo$TumorID[i] <- readFrame_clust$TumorID[1]
			clustInfo$genotype[i] <- readFrame_clust$genotype[1]
			clustInfo$library[i] <- readFrame_clust$library[1]
			clustInfo$strand[i] <- readFrame_clust$ori[1]
		}
	}

	return(clustInfo)
}	

importReadFrame <- function (dataFile, filterULflag=F) {
	
	# Requirements (for now) for determining end coordinates of chromosomes	
	library(KCsmart)
	data(mmMirrorLocs)
	
	# Cumulative chromosome end coordinates
	chromEnds <- cumsum(unlist(lapply(mmMirrorLocs, max)))
	chromEnds <- c(0, chromEnds)
	
	readFrame <- read.delim(dataFile, stringsAsFactors=F, header=T)
	
  filterInfo <- list(rawReads = nrow(readFrame))
	# Filtering

	# Filter out inserts longer than 400 bp
  longIns <- which(readFrame$end - readFrame$start > 400)
  filterInfo$longIns <- length(longIns)
	readFrame <- readFrame[-longIns,]
	# Remove chrY insertions
  chrYreads <- grep('chrY', readFrame$chr)
  filterInfo$chrYreads <- length(chrYreads)
  
	if (length(chrYreads) > 0) {
		readFrame <- readFrame[-which(readFrame$chr == 'chrY'),]
	}	
	# Remove NA barcodes
  NAbc <- grep('NA', readFrame$index)
  readFrame <- readFrame[-NAbc,]
  filterInfo$NAbc <- length(NAbc)
	# Extract a numeric version of the barcode
	readFrame$index_numeric <- as.numeric(gsub('.+BC_', '', readFrame$index))
	# Remove shearing reads that have no T7 sequence, and therefore the 
	# ligation point is unknown
  if (filterULflag) {
    noUL <- which(readFrame$flag != 1)
    readFrame <- readFrame[-noUL,]
    filterInfo$noUL <- length(noUL)
  }
  else {
    filterInfo$noUL <- 0
  }
	
	
	# Ordering and adding information
	
	readFrame$datafile <- dataFile
	# Order the reads on barcode, chromosome and startposition
	readFrame <- readFrame[order(readFrame$index_numeric, readFrame$chr, readFrame$start ),]
	# Get a numeric version of the chromosome column
	readFrame$numchr <- as.numeric(gsub('chr', '', gsub('X', '20', readFrame$chr)))
	# Assign genome-wide coordinates to insertions
	readFrame$start_con <- readFrame$start + chromEnds[readFrame$numchr]
	readFrame$end_con <- readFrame$end + chromEnds[readFrame$numchr]
	
  # Passing the filterInfo as an attribute
  
  attr(readFrame, 'filterInfo') <- filterInfo
  
	return(readFrame)
	
}	

readFrameTumorAnnotation <- function (readFrame, tumorInfo='data/tumorinfo.txt') {

	# Function to annotate a readFrame with original tumorIDs, these are
	# found in the tumorinfo.txt file for the MMTV data (which is the 
	# standard dataset)

	tumorFrame <- read.delim(tumorInfo, stringsAsFactors=F)

	uniqLib <- unique(tumorFrame$library)
	uniqDatafile <- unique(readFrame$datafile)

	if (length(uniqDatafile) > 1) {
		stop("Please only use readFrames that were constructed from a single library file")
	}	

	libName <- names(unlist(sapply(uniqLib, function(x) {grep(x, uniqDatafile)})))
  
  # Prevent occurrences of shorter library names in the long library name
  # eg: Ptenset1 occurs in Ptenset10, but we only want the longest string
  # to work with
  
  if (length(libName) > 1) {
    libName <- libName[which(nchar(libName) == max(nchar(libName)))]
  }
  
	tempTumorFrame <- tumorFrame[tumorFrame$library %in% libName,]
	tempReadFrame <- merge(x=readFrame, y=tempTumorFrame, by.x='index_numeric', by.y='BC')

  # Copy and retain attributes of the original readFrame
  filterInfo <- attr(readFrame, 'filterInfo')
  filterInfo$BCnotInTumList <- nrow(readFrame) - nrow(tempReadFrame)
  attr(tempReadFrame, 'filterInfo') <- filterInfo
  
	return(tempReadFrame)
	
}

readFrameResolveContaminations <- function (readFrame) {

	# Function to resolve contaminations in sequenced libraries. 
	# Only use on readFrames that actually have been clustered
	
  # To Do: somehow filter two highly clonal insertions as indiviual
  # insertions instead of removing them
  
  # Before contamination resolving
  noReads <- nrow(readFrame)
  
	# Define unique clusters
	uniqClust <- unique(readFrame$cluster)	

	# Check for contaminations
	contClust <- vector(mode='numeric', length=0)

	for (i in 1:length(uniqClust)) {
		readFrame_clust <- readFrame[readFrame$cluster == uniqClust[i],]
		bc <- unique(readFrame_clust$index_numeric)
		if (length(bc) > 1) {
			cat('found contamination\n')
			readFrame_cont <- readFrame[readFrame$cluster == uniqClust[i],]
			cat(paste('Cluster: ', uniqClust[i], ' - iteration ', i, '\n', sep=''))
			cat(paste('Number of unique tumors: ', length(unique(readFrame_cont$index_numeric)), '\n', sep=''))
			contClust[length(contClust)+1] <- uniqClust[i]
		}	
	}
	
	# If there are no contaminations, return the readFrame unaltered
	
	if (length(contClust) < 1) {
		cat("No contaminations found, returning readFrame\n")
		return(readFrame)
	}	
	
	# Resolve contaminations

	for (i in 1:length(contClust)) {
		readFrame_cont <- readFrame[readFrame$cluster == contClust[i],]
		# Discard all reads associated with cluster if more that 2 unique tumors
		# are associated with a cluster
    tumorDepth <- unlist(lapply(split(readFrame_cont$depth, readFrame_cont$index_numeric), sum))
		# Keep only reads from a tumor if it constitutes 80% of the reads
    # Remove the rest
    removeindex <- 
       as.numeric(names(tumorDepth[tumorDepth/sum(tumorDepth) < .8]))
		readFrame <- readFrame[-which(readFrame$cluster == contClust[i] & 
      readFrame$index_numeric %in% removeindex),]	
	}
	
  # Update filterInfo
  
  filterInfo <- attr(readFrame, 'filterInfo')
  filterInfo$contaminated <- noReads - nrow(readFrame) 
  attr(readFrame, 'filterInfo') <- filterInfo
  
  
	return(readFrame)
	
}

clustInfoCorrectClonality <- function (clustInfo) {

	# Function to correct the uniqLigation score and depth for the number
	# of reads and the total number of insertions
	
	clustInfo$corrDepth <- vector(mode='numeric', length=nrow(clustInfo))
	clustInfo$corrUniqLig <- vector(mode='numeric', length=nrow(clustInfo))
	
	uniqTum <- unique(clustInfo$TumorID)
	
	for (t in 1:length(uniqTum)) {
		totalDepth <- sum(clustInfo$depth[clustInfo$TumorID == uniqTum[t]])
		totalIns <- nrow(clustInfo[clustInfo$TumorID == uniqTum[t] & clustInfo$uniqLig > 1,])
		corrFactor <- totalDepth * totalIns
		clustInfo$corrDepth[clustInfo$TumorID == uniqTum[t]] <- clustInfo$depth[clustInfo$TumorID == uniqTum[t]] / corrFactor
		clustInfo$corrUniqLig[clustInfo$TumorID == uniqTum[t]] <- clustInfo$uniqLig[clustInfo$TumorID == uniqTum[t]] / corrFactor
	}
	
	return(clustInfo)
	
}

exportListClustInfo <- function (listClustInfo, saveFileName="insInfo.csv") {

cat(colnames(listClustInfo[[1]]), sep=',', file=saveFileName)
cat('\n',file=saveFileName, append=TRUE)

lapply(listClustInfo, write.table, quote = FALSE, sep = ",", row.names = FALSE,
	col.names=FALSE, append=TRUE, file=saveFileName)

return(1)
}
