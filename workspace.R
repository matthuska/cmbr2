# Shows the sizes of the variables in the current workspace
object.sizes <- function() {
		  return(rev(sort(sapply(ls(envir=.GlobalEnv), function (object.name) 
														  object.size(get(object.name))))))
 }
