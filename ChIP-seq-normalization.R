#' Normalizes a list of lists of numeric or data frame according a provided scheme.
#' Currently only median normalization is implemented.
#' This method uses the multicore package to parallelize computation.
#'
#' @param counts A list of lists of class numeric or data.frame.
#' @param control.vector A list of character with the same length as counts. It
#'   specifies for each list in counts, the element  from count itsself to use to
#'   normalize this element or "median"/"mean"/other function to do normalization
#'   based on distribution of the entity (e.g. DNaseI-seq) or NA if this element 
#'   should not be considered for normalization.  
#' @param method optional: Character specifying what entity normalization should 
#'   be used. Currently only "diff"(Default), "median" and "median2" are 
#'   implemented.
#' @param mc.cores Numeric indicating the number of processors to use.
#'
#' @example
#' 
#' m1 <- sample(rnorm(100, 1, .5), 50)
#' m2 <- sample(rnorm(100, 4, .2), 50)
#' m3 <- sample(rnorm(100, 4, .2), 50)
#' 
#' counts <- list("m1"=m1,"m2"=m2, "m3"=m3)
#' control.vector <- c("median", "m3")
#' normalized <- normalizeList( counts, control.vector )
#' 
#' m1 <- matrix(rnorm(100, 1, .5), ncol=20, nrow=30)
#' m2 <- matrix(rnorm(100, 4, .2), ncol=20, nrow=30)
#' m3 <- matrix(rnorm(100, 20, 3), ncol=20, nrow=30)
#' m4 <- matrix(rnorm(100, 1, .8), ncol=20, nrow=30)
#' counts <- list("m1"=m1,"m2"=m2, "m3"=m3, "m4"=m4)
#' control.vector <- c("median", "m4", "m4", NA)
#' normalized <- normalizeList( counts, control.vector )
#' 
#'
#' 
#' helmuth <helmuth@molgen.mpg.de> 2013-03-06
normalizeList <- function( counts, control.vector, method="diff", mc.cores=4 ) {
  require( multicore )
  stopifnot( method %in% c("diff", "median", "median2") )
  normalized = mclapply( 1:length(control.vector), function( i ) {
			if (is.na( control.vector[ i ]  )) { # not included (won't be an element in the output)
			  ( NA )
			} else {
			  signal = counts[[ i ]]
			  if ( is.element( control.vector[ i ] , names( counts ) )) {
			    cl= counts[[ control.vector[ i ] ]]

			    if (method == "diff" ) {
			      #helmuth 2013-11-14 RPKM_mod / RPKM_control
			      x = signal / (cl + 1)
			    } else if (method == "median" ) {

			      #helmuth for RPKM not applicable
			      r = (cl + 1) * median( (signal + 1) /  (cl + 1 ))
			      x = (signal+1) / r 

			    } else if (method == "median2") {
			      # helmuth 2013-09-04: Alternative method using less pseudo counts
			      r = cl * median( (signal + 1) /  (cl + 1 ))
			      x = signal / r 
			      x[ which(is.na(x) | is.infinite(x)) ] = 0
			    }
			    (x)
			  } else {
			    ( signal / ( get( control.vector[ i ] )( signal ) ) )
			  }
			}
}, mc.cores=mc.cores)
  normalized          = normalized[  which( !is.na( control.vector ) ) ]
  names( normalized ) = names( counts[ which( !is.na( control.vector ) ) ] )

  return ( normalized )
}
