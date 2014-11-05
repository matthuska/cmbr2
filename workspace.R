# Shows the sizes of the variables in the current workspace
object.sizes <- function() {
		  return(rev(sort(sapply(ls(envir=.GlobalEnv), function (object.name) 
														  object.size(get(object.name))))))
}

# remove outliers from a vector of numerics - "useful for plotting" (says Ho-Ryun)
remove_outliers <- function(x, na.rm = TRUE, probs=c(.05, .95), ...) {
	stopifnot(class(x) == "numeric")

  qnt <- quantile(x, na.rm = na.rm, probs=probs, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

notify <- function(title="R script info", text="Script has finished with zero exit status") {
  #system( "aplay /usr/share/sounds/pop.wav" )
  system( paste("zenity --title=\"", title, "\" --text=\"", text, "\" --info", sep="") )
}

# MAP2UNIQUE function
#returns a list with two elements:
#1: values, sorted unique values of the counts vector
#2: map, reference to the values vector such that:
# for vectors: all(values[map] == counts) 
# for matrices: all(values[map,] == counts) 
#
# Select a subsect of values: all(values[map[selection]] == counts[selection])
map2unique <- function(counts){
	if (is.matrix(counts)) {
		#Run length encoding for a matrix (groups identical rows)
		rle.matrix <- function(x) {
			n <- dim(x)[1]
			if (n == 0L) return( list(values = numeric(), map = integer()) )
			y <- !sapply(2:n, function(i) { all(x[i,] == x[(i-1),]) } )
			i <- c(which(y | is.na(y)), n)
			structure(list(lengths = diff(c(0L, i)), values = x[i,]))
		}
		o <- do.call(order, lapply(1:NCOL(counts), function(i) counts[, i]))
		uval <- rle.matrix( counts[o,] )
		n <- dim(counts)[1]
		m <- dim(uval$values)[1]
	} else if (is.vector(counts)) {
		o <- order(counts)
		uval <- rle(counts[o])
		n <- length(counts)
		m <- length(uval$values)
	}
	values <- uval$values
	uval$values <- 1:m
	map <- integer(n)
	map[o] <- inverse.rle( uval )
	return(list(values=values, map=map))
}
