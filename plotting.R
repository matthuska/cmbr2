#' Sets par to its default
#'
#' @example
#'   par(resetPar())
resetPar <- function() {
	dev.new()
	op <- par(no.readonly = TRUE)
	dev.off()
	op
}

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
plotSignature <- function(coverage, groupings=NA, grouping.order=NA, classes=NA, classes.color=NA, property=NA, property.lab="", ordering=NA, color=gray.colors(25), scale.=T, title=NA, barpl=T, markings=NULL) {
	require(gtools)

	stopifnot(length(coverage) >= 1)

	# dimensions
	k       <- length( coverage )      # number of signature components
	k.names <- names( coverage )				# names of signature components to be plotted
	n       <- dim( coverage[[1]] )[1] # number of entities in coverage
	m       <- dim( coverage[[1]] )[2] # length of the coverage profile
	if (invalid(groupings)) {
		groupings <- as.factor(rep("", n))
	} 
	if (invalid(grouping.order)) {
		grouping.order  <- levels( groupings )				# group labels
	}
	g       <- length( grouping.order )				# number of groups
	if (! invalid(classes) && invalid(classes.color)) {
		classes.color <- rainbow(length(levels(classes)))
		names(classes.color) <- levels(classes)
	}
	if (invalid(property)){
		property <- rep(0,n)
		barpl    <- F
	}
	if (invalid(ordering) ) {
		ordering <- order( property )
	}

	# data preparation
	for (i in 1:k) { rownames(coverage[[i]]) <- c() } #no duplicated rownames
	cov.mat <- data.frame( coverage )
	if (scale.) {
		cov.mat <- scale( cov.mat, 
						 scale = unlist(lapply(1:k, function(i) { rep( median(sapply( cov.mat[,((i-1)*(m+1)):(i*m)], sd)), m) })),
						 center= unlist(lapply(1:k, function(i) { rep( mean( colMeans( cov.mat[,((i-1)*(m+1)):(i*m)] )), m) }))
						 )
	} 
	limits      <- c( min( cov.mat ), max( cov.mat ) )

	# plotting layout
	cols <- k # signatures plus one column for colorcode
	widths <- c(rep(3, k), 2)
	if (! invalid(classes) ) { cols <- cols + 1; widths <- c(.2, widths) }
	if (barpl) { cols <- cols + 1; widths <- c(3, widths)	} 
	plot.matrix <- t(sapply( (0:(g-1))*cols, function( i ) { 
							c( (i+1):(i+cols) )
}))
	plot.matrix <- cbind( plot.matrix, rep( g * cols + 1, g) ) # rightmost column will be legend
	heights <- c( rep(0, g) )
	for ( i in 1:g) {
		heights[i] <- length( which(groupings == grouping.order[i]) ) / n
	}
	layout( plot.matrix, heights=heights, widths=widths) # width of columns can be specified here

	# start drawing
	par.config <- par()
	space <- par.config$cex.main
	first.row <- T
	for ( group in grouping.order ) {

		group.order <- ordering[ which( groupings == group ) ]

		if ( barpl ) {
			if (first.row) {
				par(mar=c(0, space, 2*space,1))
			} else {
				par(mar=c(0, 1.1*space, space, 1))
			}
			barplot( property[ group.order ],  ylab=paste( length(group.order), group ), cex.lab=.9*space, xlab="", horiz=T, xaxt="n", yaxt="n", yaxs="i", xlim=range(pretty(c(0, max( property )))))
			if ( first.row ) {
				axis( 3, cex.axis=.6*space )
				title( main=property.lab, cex.main=space )
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
				abline(v=markings[j], lwd=space, lty=3, col=names(markings)[j])
			}
			if (first.row) {
				par(mar=c(1, 1, 2*space, 0))
			} else {
				par(mar=c(1, 1, space, 0))
			}
		}

		first.row <- F
	}

	# add legend
	leg <- sapply(1:length(color), function(i) { "" } )
	leg[1] <- format(limits[1], digits=2)
	leg[ length(color) / 2 ] <- format(limits[2] / 2, digits=2)
	leg[ length(color) ] <- format(limits[2], digits=2)
	par(mar=c(0,0,4*space,0))
	plot(1, type = "n", axes=FALSE, xlab="", ylab="") #pseudoplot to draw legend
	legend("top", legend=rev(leg), col=rev(color), lwd=space, border="white", bg="white", bty="o", y.i=.1, horiz=F, cex=.8*space)

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
#' @param x.bins optional: A list of x values for region_length 
#' @param groupings optional: A N length vector of factors specifying the class each row in the matrix belongs to.
#' @param grouping.order optional: A N length vector of factors specifying the class each row in the matrix belongs to.
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
plotProfile <- function( coverage, x.bins=NA, groupings=NA, grouping.order=NA, color=NA, method="median", scale.=T, mark.quartiles=T, markings=NULL, ylim=NULL,  ... ) {

	stopifnot(length(coverage) >= 1)

	# dimensions
	k       <- length( coverage )      # number of signature components
	k.names <- names( coverage )				# names of signature components to be plotted
	n       <- dim( coverage[[1]] )[1] # number of entities in coverage
	m       <- dim( coverage[[1]] )[2] # length of the coverage profile
	if (invalid(x.bins))
		x.bins <- seq(0, m-1, 1)
	if (invalid(groupings)) 
		groupings <- as.factor(rep("", n))
	if (invalid(grouping.order)) 
		grouping.order  <- levels( groupings )				# group labels
	g       <- length( grouping.order )				# number of groups
	FUN     <- get( method )
	if (invalid(color)) 
		color <- rainbow( g )

	# data preparation
	for (i in 1:k) 
		rownames(coverage[[i]]) <- c() #no duplicated rownames
	cov.mat <- data.frame( coverage )
	if (scale.) 
		cov.mat <- scale( cov.mat, 
						 scale=unlist(lapply(1:k, function(i) { rep( median(sapply( cov.mat[,((i-1)*(m+1)):(i*m)], sd)), m) })),
						 center=F)
	#center= unlist(lapply(1:k, function(i) { rep( mean( colMeans( cov.mat[,((i-1)*(m+1)):(i*m)] )), m) }))
	#)

	# get ylimits
	if ( invalid(ylim) ) {
		y.min <- Inf
		y.max <- -Inf
		for ( i in 1:k ) {
			y.values <- list()
			lower.quantiles <- list()
			higher.quantiles <- list()
			for (j in 1:g) {
				what <- which( groupings == grouping.order[j] )
				y    <- cov.mat[ what, ((i-1)*m+1):(i*m) ]
				y.values[[j]] <- apply( y, 2, FUN)

				if (mark.quartiles) {
					lower.quantiles[[j]]  <- apply( y, 2, quantile, probs=.25)
					higher.quantiles[[j]] <- apply( y, 2, quantile, probs=.75)
				}

				y.min = min( y.min, unlist( y.values), unlist(lower.quantiles))
				y.max = max( y.max, unlist(y.values), unlist(higher.quantiles))	
			}
		}
		#ylim = c(y.min, y.max)
		ylim = c(0, y.max)
	} 

	# start plotting
	for ( i in 1:k ) {
		# 1st get all data to find limits
		x.values <- x.bins
		y.values <- list()
		lower.quantiles <- list()
		higher.quantiles <- list()
		for (j in 1:g) {
			what <- which( groupings == grouping.order[j] )
			y    <- cov.mat[ what, ((i-1)*m+1):(i*m) ]
			y.values[[j]] <- apply( y, 2, FUN)

			if (mark.quartiles) {
				lower.quantiles[[j]]  <- apply( y, 2, quantile, probs=.25)
				higher.quantiles[[j]] <- apply( y, 2, quantile, probs=.75)
			}
		}

		# start drawing
		plot(1, type="n", main=k.names[i], xlim=c( min(x.values) ,max(x.values)), xaxs="i", xaxt="n",  ylim=ylim, yaxs="i",  ... )
		axis(1, at=x.values, labels=T)
		grid(nx=(length(x.values)-1))

		# add legend
		legend("topright", legend=grouping.order, fill=color, border="white", bg="white", horiz=F)

		for (j in 1:g) {
			lines( x.values , y.values[[j]] , col=color[j], lty=1, ... )

			if (mark.quartiles) {
				lines( x.values , lower.quantiles[[j]], col=color[j], lty=3)
				lines( x.values , higher.quantiles[[j]], col=color[j], lty=3)
				transparent.color <- rgb(col2rgb( color[j] )[1,1]/256, col2rgb( color[j] )[2,1]/256, col2rgb( color[j] )[3,1]/256, alpha=.2)
				polygon( c(x.values, rev(x.values)), c(lower.quantiles[[j]], rev(higher.quantiles[[j]])), col=transparent.color, border=NA)
			}

			for (j in 1:length(markings)) {
				abline(v=markings[j], lty=5, col=names(markings)[j])
			}

		}

	}


}

#' A smoothScatter that is drawn to the borders of the plot (no whitespace anymore). It's suggested to use smkey instead!
#'
#' @author your highness <helmuth@molgen.mpg.de> 2013-03-14
smoothScatter.2 <- function( x, y=NULL, colramp=colorRampPalette(c("white", blues9)), xlim=NULL, ylim=NULL, postPlotHook=NULL, ... ) {

	# create a dummy plot and color it with the color for lowest value
	xy <- xy.coords(x, y)
	if ( is.null(xlim) ) {
		xlim <- c( min(xy$x), max(xy$x))
	} 
	if ( is.null(ylim) ) {
		ylim <- c( min(xy$y),  max(xy$y))
	}
	plot(0,1, xlim=xlim, ylim=ylim, type="n", ...)
	polygon( x=c( 2*xlim[1], 2*xlim[2], 2*xlim[2], 2*xlim[1]), y=c( 2*ylim[1], 2*ylim[1], 2*ylim[2], 2*ylim[2] ), col=colramp(256)[1])

	# now add smoothscatter to this plot
	smoothScatter(x=x,y=y, colramp=colramp, add=T, postPlotHook=postPlotHook, ...)
}

#' @author Ruping Sun <ruping@molgen.mpg.de>
.smoothScatterCalcDensity1 <- function(x, nbin, bandwidth, range.x) {

	if (length(nbin) == 1)
		nbin <- c(nbin, nbin)
	if (!is.numeric(nbin) || (length(nbin)!=2))
		stop("'nbin' must be numeric of length 1 or 2")

	if (missing(bandwidth)) {
		bandwidth <- diff(apply(x, 2, quantile, probs=c(0.05, 0.95), na.rm=TRUE)) / 25
	} else {
		if(!is.numeric(bandwidth))
			stop("'bandwidth' must be numeric")
	}
	## create density map
	if(missing(range.x))
		rv <- bkde2D(x, gridsize=nbin, bandwidth=bandwidth)
	else
		rv <- bkde2D(x, gridsize=nbin, bandwidth=bandwidth, range.x=range.x) 
	rv$bandwidth <- bandwidth
	return(rv)
}

#' Legend can now be suppressed with legend argument
#' 
#' @author Ruping Sun <ruping@molgen.mpg.de>
#' @author Johannes Helmuth <helmuth@molgen.mpg.de>
#'
image.plot2 = function (..., add = FALSE, nlevel = 64, legend=T, legend.shrink = 0.9, 
						legend.width = 1.2, legend.mar = NULL, graphics.reset = FALSE, 
						horizontal = FALSE, bigplot = NULL, smallplot = NULL, legend.only = FALSE, 
						col = tim.colors(nlevel)) 
{
	old.par <- par(no.readonly = TRUE)
	info <- image.plot.info(...)
	if (add) {
		big.plot <- old.par$plt
	}
	if (legend &legend.only) {
		graphics.reset <- TRUE
	}
	if (legend & is.null(legend.mar)) {
		legend.mar <- ifelse(horizontal, 3.1, 5.1)
	}
	temp <- image.plot.plt(add = add, legend.shrink = legend.shrink, 
						   legend.width = legend.width, legend.mar = legend.mar, 
						   horizontal = horizontal, bigplot = bigplot, smallplot = smallplot)
	smallplot <- temp$smallplot
	bigplot <- temp$bigplot
	if (!legend.only) {
		if (legend & !add) {
			par(plt = bigplot)
		}
		##par(bty = 'n')
		image(..., add = add, col = col)
		box()

		big.par <- par(no.readonly = TRUE)

	}
	if (legend) {
		if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
			par(old.par)
			stop("plot region too small to add legend\n")
		}
		ix <- 1
		minz <- info$zlim[1]
		maxz <- info$zlim[2]
		binwidth <- (maxz - minz)/nlevel
		midpoints <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
		iy <- midpoints
		iz <- matrix(iy, nrow = 1, ncol = length(iy))
		breaks <- list(...)$breaks
		par(new = TRUE, pty = "m", plt = smallplot, err = -1)
		if (!horizontal) {
			if (is.null(breaks)) {
				image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
					  ylab = "", col = col)
				axis(4, mgp = c(3, 1, 0), las = 2)
			}
			else {
				image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
					  ylab = "", col = col, breaks = breaks)
				axis(4, at = breaks, labels = format(breaks), mgp = c(3, 
																	  1, 0), las = 2)
			}
		}
		else {
			if (is.null(breaks)) {
				image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
					  ylab = "", col = col)
				axis(1, mgp = c(3, 1, 0))
			}
			else {
				image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
					  ylab = "", col = col, breaks = breaks)
				axis(1, at = breaks, labels = format(breaks), mgp = c(3, 
																	  1, 0))
			}
		}
		box()
		mfg.save <- par()$mfg
		if (graphics.reset | add) {
			par(old.par)
			par(mfg = mfg.save, new = FALSE)
			invisible()
		}
		else {
			par(big.par)
			par(plt = big.par$plt, xpd = TRUE)
			par(mfg = mfg.save, new = FALSE)
			invisible()
		}
	}
}

#' SmoothScatter with color key
#'
#' Legend can now be suppressed with legend argument
#'
#' @value returns the density of the ScatterPlot
#' 
#' @author Ruping Sun <ruping@molgen.mpg.de>
#' @author Johannes Helmuth <helmuth@molgen.mpg.de>
#'
smkey <- function(x, y=NULL, 
				  nbin=128,
				  bandwidth,
				  colramp=colorRampPalette(c("white", brewer.pal(9, "Blues"))),
				  legend=T, 
				  nrpoints=100,
				  transformation=function(x) x^.25,
				  xlab=NULL, ylab=NULL, postPlotHook=box,
				  pch=".", cex=1,
				  xlim, ylim, col="black",
				  xaxs=par("xaxs"), yaxs=par("yaxs"), ...) {

	if (!is.numeric(nrpoints) | (nrpoints<0) | (length(nrpoints)!=1) )
		stop("'nrpoints' should be numeric scalar with value >= 0.")

	## similar as in plot.default
	xlabel <- if (!missing(x)) 
		deparse(substitute(x))
	ylabel <- if (!missing(y)) 
		deparse(substitute(y))
	xy <- xy.coords(x, y, xlabel, ylabel)
	xlab <- if (is.null(xlab)) 
		xy$xlab
	else xlab
	ylab <- if (is.null(ylab)) 
		xy$ylab
	else ylab


	## eliminate NA
	x <- cbind(xy$x, xy$y)[!(is.na(xy$x)|is.na(xy$y)), ]

	## xlim and ylim
	if(!missing(xlim)) {
		stopifnot(is.numeric(xlim), length(xlim)==2, !any(is.na(xlim)))
		x <- x[ (x[,1]>=xlim[1]) & (x[,1]<=xlim[2]), ]
	} else {
		xlim <- range(x[,1], na.rm=TRUE)
	}
	if(!missing(ylim)) {
		stopifnot(is.numeric(ylim), length(ylim)==2, !any(is.na(ylim)))
		x <- x[ (x[,2]>=ylim[1]) & (x[,2]<=ylim[2]), ]
	} else {
		ylim <- range(x[,2], na.rm=TRUE)
	}

	## create density map
	map  <- .smoothScatterCalcDensity1(x, nbin, bandwidth, list(xlim, ylim))
	xm   <- map$x1
	ym   <- map$x2
	dens <- map$fhat
	dens <- array(transformation(dens), dim=dim(dens))	

	## plot color image
	image.plot2(xm, ym, z=dens, legend=legend, legend.shrink = 1.0,
				xlab = xlab, ylab = ylab, nlevel = 256, xlim=xlim, ylim=ylim, ...)
	if(!is.null(postPlotHook)) postPlotHook()

	# plot selection of dots
	if (nrpoints!=0){
		## we assume that map$x1 and map$x2 go linearly from
		## their first to their last value in nbin steps
		stopifnot(length(xm)==nrow(dens), length(ym)==ncol(dens))
		ixm <- round((x[,1]-xm[1])/(xm[length(xm)]-xm[1])*(length(xm)-1))
		iym <- round((x[,2]-ym[1])/(ym[length(ym)]-ym[1])*(length(ym)-1))
		idens <- dens[1 + iym*length(xm) + ixm]
		nrpoints <- min(nrow(x), ceiling(nrpoints))
		sel <- order(idens, decreasing=FALSE)[1:nrpoints]
		points(x[sel,1:2], pch=pch, cex=cex, col=col)
	}

	invisible(dens)
}

scatter.classes <- function( x, y, main="", xlab=expression(paste(plain(log)[10], "( ", plain(RPKM), " )")),
							ylab=expression(paste(plain(log)[10], "( ", plain(RPKM)["CAGE"], " )")), classes=T, nbins=1024, border=NA ) {
	source("computation/3.promoter_histone_modification_patterns/plotting.R")

	# Remove upper quartiles from x and y
	selected = NULL
	for ( clazz in c("uHCP", "dHCP", "uLCP", "dLCP")) {
		selected = c(selected, which( x[ gr$promoter.class == clazz ] < quantile(x[ gr$promoter.class == clazz ],.95) & y[ gr$promoter.class == clazz ] < quantile(y[ gr$promoter.class == clazz ],.95)))
	}
	x = x[selected]
	y = y[selected]
	granges = gr[selected]

	min.x = min( x )
	max.x = max( x )

	min.y = min( y )
	max.y = max( y )


	zlim = c(Inf, -Inf)
	for ( clazz in c("uHCP", "dHCP", "uLCP", "dLCP")) {
		map  <- .smoothScatterCalcDensity1( cbind(x[ granges$promoter.class == clazz ], y[granges$promoter.class == clazz]), nbin=nbins)
		dens <- map$fhat
		dens <- array(dens^.25, dim=dim(dens))	

		if ( min(dens) < zlim[1]) zlim[1] = min(dens)
		if ( max(dens) > zlim[2]) zlim[2] = max(dens)
	}

	layout( matrix(c(1, 2, 7, 1, 3, 7, 1, 4,  7, 5, 6, 7), ncol=3, nrow=4, byrow=T), widths=c(3,1.2,.4), height=c(3,3,3,1.2))

	#if (is.na(border)){
	#  border = max( x[ values(tss)$expression.class == "inactive"] )
	#}
	for ( clazz in c("uHCP", "dHCP", "uLCP", "dLCP")) {
		r = cor( x[ granges$promoter.class == clazz ], y[ granges$promoter.class == clazz ], method="spearman")
		print( paste(clazz, "s | ", main, " | Correlation = ", r, sep=""))
		par(mar=c(0,4,3,0))
		trash = smkey( x[ granges$promoter.class == clazz ], y[ granges$promoter.class == clazz ], main=paste(clazz, "s | SCC ", format(r, digits=2), sep=""), xlab="", ylab="", xlim=c(min.x, max.x), ylim=c(min.y, max.y), legend=F, zlim=zlim, xaxt="n")
		grid()
		if (classes)
			abline(v=border, col="white", lty=2)
		mtext(ylab, side=2, line=2)
	}

	# add title
	title(paste(main), outer=T)

	# print color bars
	par(mar=c(5,4,0,0))
	plot(1, xlim=c(min.x, max.x), ylim=c(0,1), type = "n", axes=FALSE, xlab="", ylab="", xaxs="i")
	if (classes) {
		rect( par()$usr[1],  par()$usr[3]+.1, border,        par()$usr[3]+.8, col=classes.color[1], border="white")
		rect( border,        par()$usr[3]+.1, par()$usr[2],  par()$usr[3]+.8, col=classes.color[2], border="white")
	}
	axis(1)
	mtext(xlab, side=1, line=2.5)

	plot(1, xlim=c(min.x, max.x), ylim=c(0,1), type = "n", axes=FALSE, xlab="", ylab="", xaxs="i")
	if (classes) {
		rect( par()$usr[1],  par()$usr[3]+.1, border,        par()$usr[3]+.8, col=classes.color[1], border="white")
		rect( border,        par()$usr[3]+.1, par()$usr[2],  par()$usr[3]+.8, col=classes.color[2], border="white")
	}
	axis(1)
	mtext(xlab, side=1, line=2.5)

	#legend
	old.par = par()	
	par( mar=c(1,1,1,5))
	plot(1, type = "n", axes=FALSE, xlab="", ylab="") #pseudoplot to draw legend
	image.plot(cbind( x, y ), z=zlim,legend.only=T, add=F, legend.width=15, legend.shrink=.7)

	par(old.par)
}

reducedSmoothScatter <- function( x, y, quant=.99, colramp=colorRampPalette(c("white", blues9)), part=.1, ...) {
	if ( length(x) != length(y) )
		return
	part = sample( 1:length(x), length(x)*part)
	x.sub = x[part]
	y.sub = y[part]
	selected = which( x.sub < quantile(x.sub, quant, na.rm=T) & y.sub < quantile(y.sub, quant, na.rm=T) )
	smoothScatter( x.sub[ selected ], y.sub[ selected ], colramp=colramp, ...) 
}

reducedLogSmoothScatter <- function( x, y, quant=.99, colramp=colorRampPalette(c("white", blues9)), part=.1, ...) {
	if ( length(x) != length(y) ) 
		return
	part = sample( 1:length(x), length(x)*part)
	x.sub = x[part]
	y.sub = y[part]
	selected = which( x.sub > 0 & y.sub > 0 )
	x.sub = log10( x.sub )
	y.sub = log10( y.sub )

	smoothScatter( x.sub[ selected ], y.sub[ selected ], colramp=colramp, xlab=deparse(substitute(x)),
				  ylab=deparse(substitute(y)), ...)
}
