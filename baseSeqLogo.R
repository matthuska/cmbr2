letterA <- function (x.pos, y.pos, ht, wt) 
{
	x <- c(0,  4,  6, 10,  8,  5, 2, 0, NA,2.2,2.6,7.4,7.8,2.2)
	y <- c(0, 10, 10,  0,  0,7.5, 0, 0, NA,  3,  4,  4,  3,  3)
	x <- 0.1 * x
	y <- 0.1 * y
	x <- x.pos + wt * x
	y <- y.pos + ht * y
	fill <- c("green4", "green4")
	list(x = x, y = y, fill = fill)
}

letterT <- function (x.pos, y.pos, ht, wt) 
{
	x <- c(0, 10, 10, 6, 6, 4, 4, 0)
	y <- c(10, 10, 9, 9, 0, 0, 9, 9)
	x <- 0.1 * x
	y <- 0.1 * y
	x <- x.pos + wt * x
	y <- y.pos + ht * y
	fill <- "red"
	list(x = x, y = y, fill = fill)
}

letterC <- function (x.pos, y.pos, ht, wt) 
{
	angle1 <- seq(0.3 + pi/2, pi, length = 100)
	angle2 <- seq(pi, 1.5 * pi, length = 100)
	x.l1 <- 0.5 + 0.5 * sin(angle1)
	y.l1 <- 0.5 + 0.5 * cos(angle1)
	x.l2 <- 0.5 + 0.5 * sin(angle2)
	y.l2 <- 0.5 + 0.5 * cos(angle2)
	x.l <- c(x.l1, x.l2)
	y.l <- c(y.l1, y.l2)
	x <- c(x.l, rev(x.l))
	y <- c(y.l, 1 - rev(y.l))
	x.i1 <- 0.5 + 0.35 * sin(angle1)
	y.i1 <- 0.5 + 0.35 * cos(angle1)
	x.i1 <- x.i1[y.i1 <= max(y.l1)]
	y.i1 <- y.i1[y.i1 <= max(y.l1)]
	y.i1[1] <- max(y.l1)
	x.i2 <- 0.5 + 0.35 * sin(angle2)
	y.i2 <- 0.5 + 0.35 * cos(angle2)
	x.i <- c(x.i1, x.i2)
	y.i <- c(y.i1, y.i2)
	x1 <- c(x.i, rev(x.i))
	y1 <- c(y.i, 1 - rev(y.i))
	x <- c(x, rev(x1))
	y <- c(y, rev(y1))
	x <- x.pos + wt * x
	y <- y.pos + ht * y
	fill <- "blue"
	list(x = x, y = y, fill = fill)
}

letterG <- function (x.pos, y.pos, ht, wt) 
{
	C <- letterC(0, 0, 1, 1)
	x <- C$x
	y <- C$y
	r1 <- max(x)
	h1 <- 0.4
	x <- c(x, NA, r1, 0.5, 0.5, r1 - 0.2, r1 - 0.2, r1, r1)
	y <- c(y, NA, h1, h1, h1 - 0.1, h1 - 0.1, 0, 0, h1)
	x <- x.pos + wt * x
	y <- y.pos + ht * y
	fill <- c("orange", "orange")
	list(x = x, y = y, fill = fill)
}

pwm2ic <- function(pwm) 
{
	npos <- ncol(pwm)
	ic <- numeric(length = npos)
	for (i in 1:npos) {
		ic[i] <- 2 + sum(sapply(pwm[, i], function(x) {
			if (x > 0) {
				x * log2(x)
			} else {
				0
			}
		}))
	}
	ic
}

addLetter <- function (letters, which, x.pos, y.pos, ht, wt) 
{
	if (which == "A") {
		letter <- letterA(x.pos, y.pos, ht, wt)
	}
	else if (which == "C") {
		letter <- letterC(x.pos, y.pos, ht, wt)
	}
	else if (which == "G") {
		letter <- letterG(x.pos, y.pos, ht, wt)
	}
	else if (which == "T") {
		letter <- letterT(x.pos, y.pos, ht, wt)
	}
	else {
		stop("which must be one of A,C,G,T")
	}
	letters$x <- c(letters$x, NA, letter$x)
	letters$y <- c(letters$y, NA, letter$y)
	letters$col <- c(letters$col, letter$fill)
	letters
}


baseSeqLogo <- function (pwm, ic.scale = TRUE, xlab = "Position", ylim = NULL, ...) 
{
	if (class(pwm) == "data.frame") {
		pwm <- as.matrix(pwm)
	}
	else if (class(pwm) != "matrix") {
		stop("pwm must be of class matrix or data.frame")
	}
	if (any(abs(1 - apply(pwm, 2, sum)) > 0.01)){ 
		warning("Columns of PWM should add up to 1.0, normalizing column counts")
		pwm <- apply(pwm, 2, function(co) co/sum(co))
	}
	chars <- c("A", "C", "G", "T")
	letters <- list(x = NULL, y = NULL, col = NULL)
	npos <- ncol(pwm)
	if (ic.scale) {
		ylab <- "bits"
		facs <- pwm2ic(pwm)
		if (is.null(ylim)) ylim <- 0.95*max(facs) + 0.05
	}
	else {
		ylim <- 1
		ylab <- "Probability"
		facs <- rep(1, npos)
	}
	wt <- 1
	x.pos <- 0.5
	for (j in 1:npos) {
		column <- pwm[, j]
		hts <- 0.95 * column * facs[j]
		letterOrder <- order(hts)
		y.pos <- 0
		for (i in 1:4) {
			letter <- chars[letterOrder[i]]
			ht <- hts[letterOrder[i]]
			if (ht > 0){ 
					letters <- addLetter(letters, letter, x.pos, y.pos, ht, wt)
				}
			y.pos <- y.pos + ht + 0.01
		}
		x.pos <- x.pos + wt
	}
	plot(NA, xlim=c(1,x.pos-0.5), ylim=c(0,ylim), xlab=xlab, ylab=ylab, frame.plot=F, ...)
	polygon(letters, col=letters$col, border=F)
}
