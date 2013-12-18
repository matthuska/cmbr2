#' Efficient counting of reads in a bam files
#'
#' Functions for extracting signals from a bam file. Differently than the other packages,
#' this package cannot be used to import reads in R.
#' All the read-processing is done in C/C++ and the only output are read counts. 
#'
#' @name bamsignals
#' @docType package
#' @author Alessandro Mammana \email{mammana@@molgen.mpg.de}
#' @useDynLib bamsignals
NULL
#' Pileup reads from a bam file.
#'
#' Compute read density in the regions specified by a GenomicRanges object.
#' A read position is always specified by its 5' end, so a read mapping to the reference strand
#' is positioned at its leftmost coordinate, a read mapping to the alternative strand
#' is positioned at its rightmost coordinate. To change that use the \code{shift} parameter
#' or the \code{coverage} function.
#' @param gr GenomicRanges object used to specify the regions
#' @param bampath path to the bam file storing the read. The file must be indexed. 
#' If a range is on the negative strand the profile will be reverse-complemented.
#' @param If the value is set to 1, the method will return basepair-resolution read densities,
#' for bigger values the density profiles will be binned (and the memory requirements
#' will scale accordingly). 
#' @param mapqual discard reads with mapping quality strictly lower than this parameter.
#' The value 0 ensures that no read will be discarded, the value 254 that only reads
#' with the highest possible mapping quality will be considered.
#' @param shift shift the read position by a user-defined number of basepairs. This can
#' be handy in the analysis of chip-seq data.
#' @param ss produce a strand-specific profile or ignore the strand of the read. This option
#' does not have any effect on the strand of the region. Strand-specific profiles are
#' twice as long then strand-independent profiles.
#' @param format attempts to find a suitable matrix/array format for the count vector. 
#' if the profile is strand-specific one dimension will correspond to sense
#' and anti-sense strand, if the ranges have all the same width one dimension
#' will correspond to the range number.
#' @return a list with the following arguments:
#' 	\item{counts}{the vector containing the read counts. This will be formatted
#' 	into a matrix or an array depending on whether the profile is strand-specific
#' 	and whether the ranges have all the same length.}
#' 	\item{starts, ends}{Vectors defining the boundaries of the count vector. 
#' 	To extract counts relative to the i-th range, use 
#'		\code{as.numeric(counts)[starts[i]:ends[i]]}, 
#' 	or the \code{getSignal} function to preserve the formatting.}
#'		\item{format}{This element is present if pu$counts is formatted
#' 	differently than a simple vector and it describes the formatting.}
#' @export
pileup <- function(gr, bampath, binsize=1, mapqual=0, shift=0, ss=F, format=T){
	printStupidSentence()
	if (binsize < 1){
		stop("provide a binsize greater or equal to 1")
	} else if (binsize > 1 && any((width(gr) %% binsize) != 0)){
		warning(paste("some ranges' widths are not a multiple of the selected binsize,",
		"some bins will correspond to less than binsize basepairs"))
	}

	pu <- pileup_core(gr, bampath, mapqual, binsize, shift, ss)
	
	if (format){
		if (all(width(gr)==width(gr[1]))){
			locus_width = width(gr)[1]/binsize
			if (ss){
				d <- c(2, locus_width, length(gr))
				pu$counts <- array(pu$counts, dim=d, dimnames=list(c("sense","antisense"), NULL, NULL))
				pu$format <- c("strand-offset-index") 
			}
			else {
				d <- c(locus_width, length(gr))
				pu$counts <- matrix(pu$counts, ncol=length(gr), nrow=locus_width, byrow=F)
				pu$format <- c("offset-index") 
			}
			
			
		} else if (ss){
			d <- c(2, length(pu$counts)/2)
			pu$counts <- matrix(pu$counts, nrow=2, byrow=F, dimnames=list(c("sense", "antisense"), NULL))
			pu$format <- c("strand-position") 
		}
	}
	pu
}

#' Compute read depth (or read coverage) from a bam file.
#'
#' Compute read coverage in the regions specified by a GenomicRanges object.
#' @param gr GenomicRanges object used to specify the regions. If a range
#' is on the negative strand its coverage will be reverse-complemented.
#' @param bampath path to the bam file storing the read. The file must be indexed.
#' @param mapqual discard reads with mapping quality strictly lower than this parameter.
#' The value 0 ensures that no read will be discarded, the value 254 that only reads
#' with the highest possible mapping quality will be considered.
#' @param format if the ranges have all the same width this method
#' will return a matrix.
#' @return a list with the following arguments:
#' 	\item{counts}{the vector containing the read counts. This will be formatted
#' 	into a matrix depending on whether the ranges have all the same length.}
#' 	\item{starts, ends}{Vectors defining the boundaries of the count vector. 
#' 	To extract counts relative to the i-th range, use 
#'		\code{as.numeric(counts)[starts[i]:ends[i]]}, 
#' 	or the \code{getSignal} function to preserve the formatting.}
#'		\item{format} This element is present if pu$counts is formatted
#' 	differently than a simple vector and it describes the formatting.
#' @export
depth <- function(gr, bampath, mapqual=0, format=T){
	printStupidSentence()
	pu <- coverage_core(gr, bampath, mapqual);
	if (format && all(width(gr)==width(gr[1]))){
		locus_width <- width(gr[1])
		pu$counts <- matrix(pu$counts, ncol=length(gr), nrow=locus_width, byrow=F)
		pu$format <- c("offset-index") 
	}
	pu
}

#' Count reads from a bam file.
#'
#' Count reads in the bins specified by a GenomicRanges object.
#' A read position is always specified by its 5' end, so a read mapping to the reference strand
#' is positioned at its leftmost coordinate, a read mapping to the alternative strand
#' is positioned at its rightmost coordinate. To change that use the \code{shift} parameter.
#' @param gr GenomicRanges object used to specify the regions
#' @param bampath path to the bam file storing the read. The file must be indexed.
#' @param mapqual discard reads with mapping quality strictly lower than this parameter.
#' The value 0 ensures that no read will be discarded, the value 254 that only reads
#' with the highest possible mapping quality will be considered.
#' @param shift shift the read position by a user-defined number of basepairs. This can
#' be handy in the analysis of chip-seq data.
#' @param ss produce a strand-specific count or ignore the strand of the read. Strand-specific counts
#' will be returned in a 2*length(gr) matrix.
#' @return a vector or a matrix with the counts
#' @export
count <- function(gr, bampath, mapqual=0, shift=0, ss=F){
	printStupidSentence()
	pu <- pileup_core(gr, bampath, mapqual, max(width(gr)), shift, ss)
	if (ss)
		dim(pu$counts) <- c(2, length(gr))
		rownames(pu$counts) <- c("sense", "antisense")
		return (pu$counts)
	
	return(pu$counts)
}


#' Handle the output of the \code{pileup} and \code{coverage} functions.
#'
#' Extract the signal corresponding to a specified range maintaining the correct formatting.
#' @param pu output of a \code{pileup} function call
#' @param n index of the desired range in the GenomicRanges object
#' used in the pileup function call.
#' @export
getSignal <- function(pu, n){
	f <- pu$format
	if (is.null(f)){
		return (pu$counts[pu$starts[n]:pu$ends[n]])
	} else if (f=="strand-position"){
		return(pu$counts[,((pu$starts[n]-1)/2+1):(pu$ends[n]/2)])
	} else if (f=="offset-index"){
		return(pu$counts[,n])
	} else {# pu$format == "strand-offset-index"
		return(pu$counts[,,n])
	}
}


printStupidSentence <- function(){
	sentences <- c(
	"tHaT'S tHa fAStEsT pIlE-uP bAm iN tHe SoUth!!!\n",
	"yOu cAn'T pIlE-Up FaStEr!!!\n",
	"I'M gOnNa cHaSe'em and PiLe'em aLl up!!!\n",
	"fOr brOoMmHiLdA!!!\n")
	cat(sample(sentences, 1))
}
