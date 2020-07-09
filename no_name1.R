T <- function(p, n) cos(n * acos(2 * p - 1))

make_c <- function(quantiles, p) {
  if (!all(diff(quantiles) > 0)) stop("'quantiles' must be in increasing order")
  if (!all(diff(p) > 0)) stop("'p' must be in increasing order")
  if (length(quantiles) != length(p)) stop("'quantiles' and 'p' must have the same size")
  if (any(is.infinite(quantiles))) {
    if (quantiles[1] > -Inf) {
      quantiles <- tanh(log(quantiles))
    } else if (quantiles[length(quantiles)] == Inf) {
      quantiles <- tanh(quantiles)
    } else stop("implement this")
  }
  return(solve(sapply(0:(length(p) - 1L), T, p = p), quantiles))
}

qno_name1 <- function(quantiles,
                      p = 0.5 + 0.5 * cos(seq(from = pi, to = 0, by = -pi / (length(quantiles) - 1L))),
                      check = TRUE) {
  
  if (any(is.infinite(quantiles))) {
    if (quantiles[1] > -Inf) {
      lower_bounded <- TRUE
      upper_bounded <- unbounded <- FALSE
      quantiles <- tanh(log(quantiles))
    } else if (quantiles[length(quantiles)] == Inf) {
      unbounded <- TRUE
      lower_bounded <- upper_bounded <- FALSE
      quantiles <- tanh(quantiles)
    } else stop("implement this")
  } else {
    lower_bounded <- upper_bounded <- unbounded <- FALSE
  }

  c <- make_c(quantiles, p)
  if (missing(p)) c <- zapsmall(c)
  len <- length(c) - 1L
  if (check) {
    root <- try(uniroot(function(p) {
      D <- rep(NA_real_, len + 1L)
      D[1] <- 0
      D[2] <- 1
      if (len > 1) for (k in 3:length(D)) 
        D[k] <- 2 * (2 * p - 1) * D[k - 1] - D[k - 2]
      D <- 2 * 0:len * D
      return(D %*% c)
    }, lower = 0, upper = 1, extendInt = "upX", tol = sqrt(.Machine$double.eps)), silent = TRUE)
    if (is.list(root) && root$f.root < sqrt(.Machine$double.eps) && 
        root$root > 0 && root$root < 1) {
      q <- Vectorize(function(p) (T(p, n = 0:(length(c) - 1L)) %*% c)[1])
      curve(q(p), from = 0, to = 1,
            n = 10001, xname = "p", ylab = expression(theta), axes = FALSE, las = 1)
      axis(1)
      abline(v = root$root, col = 2, lty = 2)
      stop("\nImplied quantile function is decreasing just before p = ", 
           round(root$root, digits = 3), ".",
           "\nTry increasing the number of quantiles and / or changing their values.")
    }
  }
  if (lower_bounded)
    return(Vectorize(function(p) exp(atanh(crossprod(T(p, n = 0:(length(c) - 1L)), c)[1]))))
  if (unbounded)
    return(Vectorize(function(p) atanh(crossprod(T(p, n = 0:(length(c) - 1L)), c)[1])))
  
  return(Vectorize(function(p) crossprod(T(p, n = 0:(length(c) - 1L)), c)[1]))
}
