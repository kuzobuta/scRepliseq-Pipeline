#' Mappability correction
#'
#' Correct a list of \code{\link{binned.data}} by mappability.
#' 
#' @param binned.data.list A \code{list} with \code{\link{binned.data}} objects or a list of filenames containing such objects.
#' @param reference A file or \code{\link{GRanges}} with aligned reads.
#' @param same.binsize If \code{TRUE} the mappability correction will only be calculated once. Set this to \code{TRUE} if all \code{\link{binned.data}} objects describe the same genome at the same binsize.
#' @return A \code{list} with \code{\link{binned.data}} objects with adjusted read counts.
#' @author Aaron Taudt
#' @inheritParams bam2GRanges
#' @inheritParams bed2GRanges
#'
#' This script was originally from "correctionMethods.R" of AneuFinder (1.2.1)
#' Aaron Taudt <aaron.taudt@gmail.com>
#' Extraction from the original script was done by Hisashi
#' 

correctMappability <- function(binned.data.list, same.binsize, reference, assembly, pairedEndReads=FALSE, min.mapq=10, remove.duplicate.reads=TRUE, max.fragment.width=1000) {
  
  binned.data.list <- loadFromFiles(binned.data.list, check.class='GRanges')
  same.binsize.calculated <- FALSE
  for (i1 in 1:length(binned.data.list)) {
    binned.data <- binned.data.list[[i1]]
    
    ## Calculate GC content per bin
    if (same.binsize & !same.binsize.calculated | !same.binsize) {
      ptm <- startTimedMessage("Calculating mappability per bin ...")
      refbin <- binReads(file=reference, assembly=assembly, chromosomes=seqlevels(binned.data), pairedEndReads=pairedEndReads, min.mapq=min.mapq, remove.duplicate.reads=remove.duplicate.reads, max.fragment.width=max.fragment.width, binsizes=NULL, reads.per.bin=NULL, bins=list('ref'=binned.data), save.as.RData=FALSE, calc.complexity=FALSE)[[1]]
      ## Check if seqlengths of data and mappability correction are consistent
      chromlengths <- seqlengths(binned.data)
      chroms <- names(chromlengths)
      # Compare
      compare <- chromlengths[chroms] == seqlengths(refbin)[chroms]
      if (any(compare==FALSE, na.rm=TRUE)) {
        warning(paste0(attr(binned.data,'ID'),": Chromosome lengths differ between binned data and 'reference'. Mappability correction skipped. Please use the correct genome for option 'reference'."))
        binned.data.list[[i1]] <- binned.data
        next
      }
      
      ## Make the mappability correction vector
      tab <- table(refbin$counts)
      refbin.maxcount <- as.numeric(names(which.max(tab[as.numeric(names(tab))>0])))
      mappability <- refbin$counts / refbin.maxcount
      mappability[mappability==0] <- 1
      
      
      same.binsize.calculated <- TRUE
      stopTimedMessage(ptm)
    }
    binned.data$mappability <- mappability
    
    ### GC correction ###
    ptm <- startTimedMessage("Mappability correction ...")
    counts <- binned.data$counts / binned.data$mappability
    mcounts <- binned.data$mcounts / binned.data$mappability
    pcounts <- binned.data$pcounts / binned.data$mappability
    ## Correction factors
    binned.data$counts <- as.integer(round(counts))
    binned.data$mcounts <- as.integer(round(mcounts))
    binned.data$pcounts <- as.integer(round(pcounts))
    binned.data$counts[binned.data$counts<0] <- 0
    binned.data$pcounts[binned.data$pcounts<0] <- 0
    binned.data$mcounts[binned.data$mcounts<0] <- 0
    
    ### Quality measures ###
    ## Spikyness
    attr(binned.data, 'spikiness') <- qc.spikiness(binned.data$counts)
    ## Shannon entropy
    attr(binned.data, 'entropy') <- qc.entropy(binned.data$counts)
    
    binned.data.list[[i1]] <- binned.data
  }
  return(binned.data.list)
}

startTimedMessage <- function(...) {
  
  x <- paste0(..., collapse='')
  message(x, appendLF=FALSE)
  ptm <- proc.time()
  return(ptm)
  
}


stopTimedMessage <- function(ptm) {
  
  time <- proc.time() - ptm
  message(" ", round(time[3],2), "s")
  
}

#' Quality control measures for binned read counts
#'
#' Calculate various quality control measures on binned read counts.
#'
#' The Shannon entropy is defined as
#' \eqn{S = - sum( n * log(n) )}, where \eqn{n = counts/sum(counts)}.\cr\cr 
#' Spikyness is defined as \eqn{K = sum(abs(diff(counts))) / sum(counts)}.
#' 
#' @param counts A vector of binned read counts.
#' @param hmm An \code{\link{aneuHMM}} object.
#' @return A numeric.
#' @name qualityControl
#' @author Aaron Taudt
#' @describeIn qualityControl Calculate the spikiness of a library

qc.spikiness <- function(counts) {
  if (is.null(counts)) {
    return(NA)
  }
  counts <- as.vector(counts)
  sum.counts <- sum(counts)
  spikiness <- sum(abs(diff(counts))) / sum.counts
  return(spikiness)
}

#' @describeIn qualityControl Calculate the Shannon entropy of a library
qc.entropy <- function(counts) {
  if (is.null(counts)) {
    return(NA)
  }
  counts <- as.vector(counts)
  total.counts <- sum(counts)
  n <- counts/total.counts
  entropy <- -sum( n * log(n) , na.rm=TRUE)
  return(entropy)
}

