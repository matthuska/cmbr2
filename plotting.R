#' Plot chromatin signature for a set of equal-sized regions
#' Plots a chromatin signature plot for a given list of equal-dimensionalized matrices with real numbers to the output device. The entities can be classified by argument and/or ordered by an argument. The color scale can be specified for convenient display. A colorcode is print on the right side of the signature plot. You can control fontsize etc. with the par() command.
#' 
#' @param coverage A list of matrices with dimensions N x region_length. The names in the list are printed above each signature plot. List can also be of length 1.
#' @param classes optional: A N length vector of factors specifying the class each row in the matrix belongs to.
#' @param ordering optional: A N length vector of ordering values. If NA, rows in matrices of coverage are plot in sequential order or as specified with property argument.
#' @param property optional: A N length vector of numeric (e.g. expression values). A barplot is drawn indicating these values for each row in coverage matrices. If ordering is NA, the property is used to sort output.
#' @param color optional: A vector of colors used for the heatmap plotting (e.g. gray.colors(25))
#' @param scale optional: If TRUE scaling and centering is applied to the data matrices in coverage. Therefore, the mean of the columnMeans of a matrix from coverage are used to center it. If FALSE, no scaling or centering is done inside the function.
#' @param title optional: A character vector printed to the outer margin of the plot as a title.
#' @param barpl optional: If TRUE a barplot with ordering values is drawn on the left side.
#' 
#' @examples
#' 
#'   m1 <- matrix(rnorm(100, 1, .5), ncol=20, nrow=30)
#'   m2 <- matrix(rnorm(100, 4, .2), ncol=20, nrow=30)
#'   m3 <- matrix(rnorm(100, 20, 3), ncol=20, nrow=30)
#'   coverage <- list("m1"=m1,"m2"=m2, "m3"=m3)
#'   classes <- as.factor( c(rep(2,10), rep(1,5), rep(3, 10), rep(2, 5)) )
#'   ordering <- 1:30
#'   property <- log( ordering )
#'   
#'   X11()
#'   plotSignature( coverage, classes, property, ordering )
#' 	
#' 
#' @author Johannes Helmuth
#' @date  2013/03/07
#' @export
plotSignature <- function(coverage, classes=NA, property=NA, ordering=NA, color=gray.colors(25), scale.=T, title=NA, barpl=T) {

	stopifnot(length(coverage) >= 1)

	# dimensions
	k       = length( coverage )      # number of signature components
	k.names = names( coverage )				# names of signature components to be plotted
	n       = dim( coverage[[1]] )[1] # number of entities in coverage
	m       = dim( coverage[[1]] )[2] # length of the coverage profile
	if (invalid(classes)) {
		classes = as.factor(rep("", n))
	} 
	groups  = levels( classes )				# group labels
	g       = length( groups )				# number of groups
	if (invalid(property)){
		property = rep(0,n)
		barpl    = F
	}
	if (invalid(ordering) ) {
		ordering = order( property )
	}

	# data preparation
	for (i in 1:k) { rownames(coverage[[i]]) = c() } #no duplicated rownames
	cov.mat = data.frame( coverage )
	if (scale.) {
		cov.mat = scale( cov.mat, 
										scale = unlist(lapply(1:k, function(i) { rep( median(sapply( cov.mat[,((i-1)*(m+1)):(i*m)], sd)), m) })),
										center= unlist(lapply(1:k, function(i) { rep( mean( colMeans( cov.mat[,((i-1)*(m+1)):(i*m)] )), m) }))
							)
	} 
	limits      = c( min( cov.mat ), max( cov.mat ) )

	# plotting layout
	cols = k # signatures plus one column for colorcode
	if (barpl) { cols = cols + 1	} 
	plot.matrix = t(sapply( (0:(g-1))*cols, function( i ) { 
										               c( (i+1):(i+cols) )
							                   }))
	plot.matrix = cbind( plot.matrix, rep( g * cols + 1, g) ) # rightmost column will be legend
	heights = c( rep(0, g) )
	for ( i in 1:length( groups )) {
			heights[i] = length( which(classes == groups[i]) ) / n
	}
	layout( plot.matrix, heights=heights) # width of columns can be specified here

	# start drawing
	par.config = par()
	space = par.config$cex.main
	first.row = T
	for ( group in groups ) {

		group.order = ordering[ which( classes == group ) ]

		if ( barpl ) {
			if (first.row) {
				par(mar=c(0, 4*space, 2*space, space))
			} else {
				par(mar=c(0, 4*space, space, space))
			}
			barplot( property[ group.order ],  ylab=paste( group ), cex.lab=space, xlab="", horiz=T, xaxt="n", yaxt="n", yaxs="i", xlim=range(pretty(c(0, max( property )))))
			if ( first.row ) {
				axis( 3 )
			}
		}

		if (first.row) {
			par(mar=c(0, space, 2*space, 0))
		} else {
			par(mar=c(0, space, space, 0))
		}
		for (i in 1:k) {
			image( t(cov.mat[ group.order, ((i-1)*m+1):(i*m) ] ), col=color, axes=F, zlim=limits)
			if ( first.row ) {
				title( main=paste( k.names[ i ] ) )  
			}
			abline(v=0.5, lwd=space, lty=3, col="black")
		}
		
		first.row = F
	}

	# add legend
	leg = sapply(1:length(color), function(i) { "" } )
	leg[1] = format(limits[1], digits=2)
	leg[ length(color) / 2 ] = format(limits[2] / 2, digits=2)
	leg[ length(color) ] = format(limits[2], digits=2)
	par(mar=c(0,2*space,2*space,0))
	plot(1, type = "n", axes=FALSE, xlab="", ylab="") #pseudoplot to draw legend
	legend("top", legend=rev(leg), col=rev(color), lwd=par()$lwd * 5, border="white", bg="white", bty="o", y.i=.1, horiz=F)

	# write title
	if (! is.na(title) ) {
		title(main=title, outer=T)
	}
}

#' Plots a coverage profile of a supplied list of coverages
plotAverageProfile <- function( coverage, classes=NA, method="mean", ... ) {
	warning("helmuth 2013-03-07: Not implemented yet!")
}


