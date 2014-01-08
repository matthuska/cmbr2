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

