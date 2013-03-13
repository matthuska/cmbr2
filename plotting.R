#' Plot chromatin signature for a set of equal-sized regions
#' Plots a chromatin signature plot for a given list of equal-dimensionalized matrices with real numbers to the output device. The entities can be classified by argument and/or ordered by an argument. The color scale can be specified for convenient display. A colorcode is print on the right side of the signature plot. You can control fontsize etc. with the par() command.
#' 
#' @param coverage A list of matrices with dimensions N x region_length. The names in the list are printed above each signature plot. List can also be of length 1.
#' @param groupings optional: A N length vector of factors specifying the class each row in the matrix belongs to. This is used to divide the signature plots in several grouping dependent signature plots.
#' @param classes optional: A N length vector of factors specifying a class each row in the matrix also belongs to. This is drawn as a color code in the left most column of the output.
#' @param classes.color optional: A named vector for the colors for the classes color code where names are classes and elements are colors. Has the same length as there are classes.
#' @param property optional: A N length vector of numeric (e.g. expression values). A barplot is drawn indicating these values for each row in coverage matrices. If ordering is NA, the property is used to sort output.
#' @param ordering optional: A N length vector of ordering values. If NA, rows in matrices of coverage are plot in sequential order or as specified with property argument.
#' @param color optional: A vector of colors used for the heatmap plotting (e.g. gray.colors(25))
#' @param scale optional: If TRUE scaling and centering is applied to the data matrices in coverage. Therefore, the mean of the columnMeans of a matrix from coverage are used to center it. If FALSE, no scaling or centering is done inside the function.
#' @param title optional: A character vector printed to the outer margin of the plot as a title.
#' @param barpl optional: If TRUE a barplot with ordering values is drawn on the left side.
#' @param markings optional: A named list of x-points where a vertical line should be drawn. The names of the list entries specify the color of the line.
#' 
#' @examples
#' 
#'   m1 <- matrix(rnorm(100, 1, .5), ncol=20, nrow=30)
#'   m2 <- matrix(rnorm(100, 4, .2), ncol=20, nrow=30)
#'   m3 <- matrix(rnorm(100, 20, 3), ncol=20, nrow=30)
#'   coverage <- list("m1"=m1,"m2"=m2, "m3"=m3)
#'   groupings <- as.factor( c(rep(2,10), rep(1,5), rep(3, 10), rep(2, 5)) )
#'   ordering <- 1:30
#'   property <- log( ordering )
#'   
#'   X11()
#'   plotSignature( coverage, groupings, property, ordering )
#' 	
#' 
#' @author Johannes Helmuth <helmuth@molgen.mpg.de>
#' @date  2013/03/07
#' @export
plotSignature <- function(coverage, groupings=NA, classes=NA, classes.color=NA, property=NA, ordering=NA, color=gray.colors(25), scale.=T, title=NA, barpl=T, markings=NULL) {

	stopifnot(length(coverage) >= 1)

	# dimensions
	k       = length( coverage )      # number of signature components
	k.names = names( coverage )				# names of signature components to be plotted
	n       = dim( coverage[[1]] )[1] # number of entities in coverage
	m       = dim( coverage[[1]] )[2] # length of the coverage profile
	if (invalid(groupings)) {
		groupings = as.factor(rep("", n))
	} 
	groups  = levels( groupings )				# group labels
	g       = length( groups )				# number of groups
	if (! invalid(classes) && invalid(classes.color)) {
		classes.color = rainbow(length(levels(classes)))
		names(classes.color) = levels(classes)
	}
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
	widths = c(rep(2, k), 1)
	if (! invalid(classes) ) { cols = cols + 1; widths = c(.2, widths) }
	if (barpl) { cols = cols + 1; widths = c(1, widths)	} 
	plot.matrix = t(sapply( (0:(g-1))*cols, function( i ) { 
										               c( (i+1):(i+cols) )
							                   }))
	plot.matrix = cbind( plot.matrix, rep( g * cols + 1, g) ) # rightmost column will be legend
	heights = c( rep(0, g) )
	for ( i in 1:length( groups )) {
			heights[i] = length( which(groupings == groups[i]) ) / n
	}
	layout( plot.matrix, heights=heights, widths=widths) # width of columns can be specified here

	# start drawing
	par.config = par()
	space = par.config$cex.main
	first.row = T
	for ( group in groups ) {

		group.order = ordering[ which( groupings == group ) ]

		if ( barpl ) {
			if (first.row) {
				par(mar=c(0, 5*space, 2*space,1))
			} else {
				par(mar=c(0, 5*space, space, 1))
			}
			barplot( property[ group.order ],  ylab=paste( length(group.order), group ), cex.lab=space, xlab="", horiz=T, xaxt="n", yaxt="n", yaxs="i", xlim=range(pretty(c(0, max( property )))))
			if ( first.row ) {
				axis( 3 )
			}
		}

		if (! invalid(classes) ) {
			if (first.row) {
				par(mar=c(0, 0, 2*space,0))
			} else {
				par(mar=c(0, 0, space, 0))
			}
			image( t(matrix(classes[ ordering ], ncol=1)), col=classes.color, axes=F)
		}

		if (first.row) {
			par(mar=c(1, 1, 2*space, 0))
		} else {
			par(mar=c(1, 1, space, 0))
		}

		for (i in 1:k) {
			image( t(cov.mat[ group.order, ((i-1)*m+1):(i*m) ] ), col=color, axes=F, zlim=limits)
			if ( i==1 & !barpl) {
				title( ylab=paste( length(group.order), group ), cex.lab=space, line=1, las=1)  
			}
			if ( first.row ) {
				title( main=paste( k.names[ i ] ) )  
			}
			for (j in 1:length(markings)) {
				abline(v=markings[j], lwd=2*space, lty=3, col=names(markings)[j])
			}
			if (first.row) {
				par(mar=c(1, 1, 2*space, 0))
			} else {
				par(mar=c(1, 1, 1, 0))
			}
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
	legend("top", legend=rev(leg), col=rev(color), lwd=space, border="white", bg="white", bty="o", y.i=.1, horiz=F, cex=space)

	# write title
	if (! is.na(title) ) {
		title(main=title, outer=T)
	}
}

#' Plots a coverage profile of a supplied list of coverages
#' 
#' Plots a chromatin signature plot for a given list of equal-dimensionalized matrices with real numbers to the output device. The entities can be classified by argument and/or ordered by an argument. The color scale can be specified for convenient display. A colorcode is print on the right side of the signature plot. You can control fontsize etc. with the par() command.
#' 
#' @param coverage A list of matrices with dimensions N x region_length. The names in the list are printed above each signature plot. List can also be of length 1.
#' @param classes optional: A N length vector of factors specifying the class each row in the matrix belongs to.
#' @param color optional: A vector of colors used for the heatmap plotting. Should be of length length( levels(classes) )
#' @param method optional: A character vector defining which method to use for condensing the coverage matrices.
#' @param scale. optional: If TRUE scaling and centering is applied to the data matrices in coverage. Therefore, the mean of the columnMeans of a matrix from coverage are used to center it. If FALSE, no scaling or centering is done inside the function.
#' @param mark.quartiles optional: If TRUE, empirical quartiles for each profile are indicated with dotted lines in the same color than the corresponding profile.
#' @param markings optional: A named list of x-points where a vertical line should be drawn. The names of the list entries specify the color of the line.
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
#'   plotProfile( coverage,   )
#' 	
#' 
#' @author Johannes Helmuth
#' @date  2013/03/11
#' @export
#' 
#' TODO helmuth 2013-03-12 provide method to set x ticks
plotProfile <- function( coverage, classes=NA, color=NA, method="median", scale.=F, mark.quartiles=T, markings=NULL, ... ) {

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
	FUN     = get( method )
	if (invalid(color)) {
		color = rainbow( g )
	}

	# data preparation
	for (i in 1:k) { rownames(coverage[[i]]) = c() } #no duplicated rownames
	cov.mat = data.frame( coverage )
	if (scale.) {
		cov.mat = scale( cov.mat, 
										scale = unlist(lapply(1:k, function(i) { rep( median(sapply( cov.mat[,((i-1)*(m+1)):(i*m)], sd)), m) }))
										#center= unlist(lapply(1:k, function(i) { rep( mean( colMeans( cov.mat[,((i-1)*(m+1)):(i*m)] )), m) }))
							)
	} 

	# start plotting
	for ( i in 1:k ) {
		# 1st get all data to find limits
		x.values = seq(1, m, 1)
		y.values = list()
		y.min = Inf
		y.max = -Inf
		lower.quantiles = list()
		higher.quantiles = list()
		for (j in 1:g) {
			what = which( classes == groups[j] )
			y    = cov.mat[ what, ((i-1)*m+1):(i*m) ]
			y.values[[j]] = apply( y, 2, FUN)


			if (mark.quartiles) {
				lower.quantiles[[j]]  = apply( y, 2, quantile, probs=c(.25, .75))[1,]
				higher.quantiles[[j]] = apply( y, 2, quantile, probs=c(.25, .75))[2,]
			}
		}

		# start drawing
		plot(1, type="n", main=k.names[i], ylim=c(0, max( c( unlist(y.values), unlist(higher.quantiles) ) )),, xlim=c(0,m), xlab="", ylab="", ... )
		for (j in 1:g) {
			lines( x.values , y.values[[j]] , col=color[j], lty=1, ... )

			if (mark.quartiles) {
				lines( x.values , lower.quantiles[[j]], col=color[j], lty=3)
				lines( x.values , higher.quantiles[[j]], col=color[j], lty=3)
				transparent.color = rgb(col2rgb( color[j] )[1,1]/256, col2rgb( color[j] )[2,1]/256, col2rgb( color[j] )[3,1]/256, alpha=.2)
				polygon( c(x.values, rev(x.values)), c(lower.quantiles[[j]], rev(higher.quantiles[[j]])), col=transparent.color, border=NA)
			}

			for (j in 1:length(markings)) {
				abline(v=markings[j], lty=5, col=names(markings)[j])
			}

			# add legend
			legend("topright", legend=groups, fill=color, border="white", bg="white", horiz=F, ...)
		}

		grid()
	}


}


