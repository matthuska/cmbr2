# author: Mike Love
# when: 2013 Feb 6
# what: some functions for speeding up calculations on sparse matrices
# as defined in the Matrix package


#' Sparse correlation, covariance, cosine similarity or Euclidean distance
#'
#' Calcultion of matrices of various similarity/distances, which takes use
#' of matrix multiplication.  Note: the comparisons are made between
#' columns, similar to the \code{cor} function, and unlike the
#' \code{dist} function which compares between rows.
#'
#' See the timing vignette for comparison of the performance of functions
#' against standard dense calculations. 
#'
#' @param X a matrix X
#' @param Y optional: a matrix in which case columns of X will be compared to
#' columns of Y.  The resulting matrix will be ncol(X) by ncol(Y).
#' @param XtX optional: if not using Y, supplying XtX, the crossproduct of X with
#' itself, will be used rather than recalculating this in the function.
#'
#' @return the correlation matrix of X, or X and Y;
#' the covariance matrix of X, or X and Y;
#' the cosine similarity matrix of X, or X and Y;
#' the Euclidean distances matrix of X, or X and Y.
#'
#' @note The sparseCot function was developed from a discussion here:
#' \url{http://stackoverflow.com/questions/5888287/running-cor-or-any-variant-over-a-sparse-matrix-in-r}
#'
#' @examples
#'
#'   sds1 <- simulateSparseData(100, c(5,5))
#'   sds2 <- simulateSparseData(100, c(2,2))
#'   cormat <- sparseCor(sparseData(sds1))
#'   cormatxy <- sparseCor(sparseData(sds1),sparseData(sds2))
#'
#' @export
sparseCor <- function(X,Y=NULL,XtX=NULL) {
  if (!is(X,"dgCMatrix")) stop("X should be a dgCMatrix")
  if (is.null(Y)) {
    if (is.null(XtX)) {
      XtX <- crossprod(X)
    } else {
      if (ncol(XtX) != ncol(X) | nrow(XtX) != ncol(X)) stop("XtX should have same number of rows and columns as number of columns of X")
    }
    n <- nrow(X)
    cMeans <- colMeans(X)
    covmat <- (as.matrix(XtX) - n*tcrossprod(cMeans))/(n-1)
    sdvec <- sqrt(diag(covmat))
    cormat <- covmat/crossprod(t(sdvec))
    return(cormat)
  } else {
    if (!is(Y,"dgCMatrix")) stop("Y should be a dgCMatrix")
    if (nrow(X) != nrow(Y)) stop("X and Y should have the same number of rows")
    n <- nrow(X)
    cMeansX <- colMeans(X)
    cMeansY <- colMeans(Y)
    covmat <- (as.matrix(crossprod(X,Y)) - n * tcrossprod(cMeansX,cMeansY))/(n-1)
    sdvecX <- sqrt(diag((as.matrix(crossprod(X)) - n*tcrossprod(cMeansX))/(n-1)))
    sdvecY <- sqrt(diag((as.matrix(crossprod(Y)) - n*tcrossprod(cMeansY))/(n-1)))
    cormat <- covmat/outer(sdvecX,sdvecY)
    return(cormat)
  }
}
sparseCov <- function(X,Y=NULL,XtX=NULL) {
  if (!is(X,"dgCMatrix")) stop("X should be a dgCMatrix")
  if (is.null(Y)) {
    if (is.null(XtX)) {
      XtX <- crossprod(X)
    } else {
      if (ncol(XtX) != ncol(X) | nrow(XtX) != ncol(X)) stop("XtX should have same number of rows and columns as number of columns of X")
    }
    n <- nrow(X)
    cMeans <- colMeans(X)
    covmat <- (as.matrix(XtX) - n*tcrossprod(cMeans))/(n-1)
    return(covmat)
  } else {
    if (!is(Y,"dgCMatrix")) stop("Y should be a dgCMatrix")
    if (nrow(X) != nrow(Y)) stop("X and Y should have the same number of rows")
    n <- nrow(X)
    cMeansX <- colMeans(X)
    cMeansY <- colMeans(Y)
    covmat <- (as.matrix(crossprod(X,Y)) - n * tcrossprod(cMeansX,cMeansY))/(n-1)
    return(covmat)
  }
}
sparseCosine <- function(X,Y=NULL,XtX=NULL){
if (!is(X,"dgCMatrix")) stop("X should be a dgCMatrix")
  if (is.null(Y)) {
    if (is.null(XtX)) {
      XtX <- crossprod(X)
    } else {
      if (ncol(XtX) != ncol(X) | nrow(XtX) != ncol(X)) stop("XtX should have same number of rows and columns as number of columns of X")
    }
    lengths <- sqrt(diag(as.matrix(XtX)))
    return(as.matrix(XtX/crossprod(t(lengths))))
  } else {
    if (!is(Y,"dgCMatrix")) stop("Y should be a dgCMatrix")
    if (nrow(X) != nrow(Y)) stop("X and Y should have the same number of rows")
    lengths.X <- sqrt(diag(as.matrix(crossprod(X))))
    lengths.Y <- sqrt(diag(as.matrix(crossprod(Y))))
    XtY <- crossprod(X,Y)
    return(as.matrix(XtY/outer(lengths.X,lengths.Y)))
  }
}
sparseEuclid <- function(X,Y=NULL,XtX=NULL) {
if (!is(X,"dgCMatrix")) stop("X should be a dgCMatrix")
  if (is.null(Y)) {
    if (is.null(XtX)) {
      XtX <- crossprod(X)
    } else {
      if (ncol(XtX) != ncol(X) | nrow(XtX) != ncol(X)) stop("XtX should have same number of rows and columns as number of columns of X")
    }
    X.length <- colSums(X^2)
    length.sum <- outer(X.length,X.length,"+")
    dimnames(length.sum) <- dimnames(XtX)
    return(as.matrix(sqrt(length.sum - 2*XtX)))
  } else {
    if (!is(Y,"dgCMatrix")) stop("Y should be a dgCMatrix")
    if (nrow(X) != nrow(Y)) stop("X and Y should have the same number of rows")
    XtY <- crossprod(X,Y)
    X.length <- colSums(X^2)
    Y.length <- colSums(Y^2)
    length.sum <- outer(X.length,Y.length,"+")
    dimnames(length.sum) <- dimnames(XtY)
    return(as.matrix(sqrt(length.sum - 2*XtY)))
  }
}

#' Apply function sparsely
#' 
#' Applies a function to the nonzero elements of a sparse vector or matrix
#' 
#' @param z a sparse vector or matrix
#' @param f a function
#'
#' @return returns the variable z, with only the value of the nonzero elements changed
#'
#' @examples
#'    sd <- simulateSparseData(100, c(5,5))
#'    logPlusOne <- function(x) log(x + 1)
#'    sparseData(sds) <- applyFunctionSparsely(sparseData(sds), logPlusOne)
#'
#' @export
applyFunctionSparsely <- function(z, f) {
  if (!(is(z,"dsparseMatrix") | is(z,"dsparseVector"))) stop("z should be a dsparseMatrix or dsparseVector")
  z@x <- f(z@x)
  z
}

#' Sparse threshold
#' 
#' Pushes the smaller values (in absolute value) to zero to achieve a
#' certain nonzero ratio.  In other words, the nonzero values closest to
#' zero are assigned zero to achieve at least the specified nonzero ratio.
#' 
#' @param z a sparse vector or matrix
#' @param nzr the desired nonzero ratio
#'
#' @return returns the variable z, with more or equal nonzero elements to achieve
#' at least a nonzero ratio of nzr
#'
#' @examples
#'   sd <- simulateSparseData(100, c(5,5), nzs=.5, nzg=.3)
#'   nnzero(sparseData(sds))/prod(dim(sds))
#'   sparseData(sds) <- sparseThreshold(sparseData(sds), .1)
#'   nnzero(sparseData(sds))/prod(dim(sds))
#'
#' @export
sparseThreshold <- function(z, nzr = .1) {
  if (!(is(z,"dsparseMatrix") | is(z,"dsparseVector"))) stop("z should be a dsparseMatrix or dsparseVector")
  if (nzr < 0 | nzr > 1) stop("nzr should be a value between 0 and 1")
  nzratio <- nnzero(z)/length(z)
  if (nzratio > nzr) {
    threshold <- quantile(abs(z@x), 1 - nzr/nzratio)
    z[z != 0 & z <= threshold & z >= -threshold] <- 0
  }
  as(z,"sparseMatrix")
}

