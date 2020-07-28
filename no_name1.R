T <- function(p, n) cos(n * acos(2 * p - 1))

make_c <- function(quantiles, u) {
  if (!all(diff(quantiles) > 0)) stop("'quantiles' must be in increasing order")
  if (!all(diff(u) > 0)) stop("'p' must be in increasing order")
  if (length(quantiles) != length(u)) stop("'quantiles' and 'p' must have the same size")
  if (any(is.infinite(quantiles))) {
    if (quantiles[1] > -Inf) {
      quantiles <- tanh(log(quantiles))
    } else if (quantiles[length(quantiles)] == Inf) {
      quantiles <- tanh(quantiles)
    } else stop("implement this")
  }
  return(solve(sapply(0:(length(u) - 1L), T, p = u), quantiles))
}

qno_name1 <- function(quantiles,
                      u = 0.5 + 0.5 * cos(seq(from = pi, to = 0, by = -pi / (length(quantiles) - 1L))),
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

  c <- make_c(quantiles, u)
  if (missing(u)) c <- zapsmall(c)
  len <- length(c) - 1L
  if (check) boyd(c)
  if (lower_bounded)
    return(Vectorize(function(p) exp(atanh(crossprod(T(p, n = 0:(length(c) - 1L)), c)[1]))))
  if (unbounded)
    return(Vectorize(function(p) atanh(crossprod(T(p, n = 0:(length(c) - 1L)), c)[1])))
  
  return(Vectorize(function(p) crossprod(T(p, n = 0:(length(c) - 1L)), c)[1]))
}

# https://epubs.siam.org/doi/pdf/10.1137/S0036142901398325 section 3
boyd <- function(a) {
  SEQ <- seq(from = 1, to = length(a), by = 2L)
  a_even <- a[ SEQ]
  a_odd  <- a[-SEQ]
  
  if (length(a_even) == 1L) Q_even <- diag(1)
  else {
    Q_even <- diag(2.0^(2L * (1L:length(a_even)) - 3L))
    Q_even[1L, 1L] <- 1.0
  }
  if (ncol(Q_even) > 1L) for (j in 2:ncol(Q_even)) {
    j2 <- 2L * j
    j4 <- 4L * j
    for (K in 1:(j - 1L)) {
      K2 <- 2L * K
      d <- j2 - K2
      Q_even[j - K, j] <- round(-d * (d - 1) /
                          (K2 * (j4 - K2 - 4)) * Q_even[j - K + 1L, j])
    }
  }
  
  Q_odd <- diag(2.0^(2L * (1L:length(a_odd)) - 2L))
  if (ncol(Q_odd) > 1L) for (j in 2:ncol(Q_odd)) {
    j2 <- 2L * j
    j4 <- 4L * j
    for (K in 1:(j - 1L)) {
      K2 <- 2L * K
      d <- j - K
      Q_odd[j - K, j] <- round(-(2 * d + 1) * d /
                              (K * (j4 - K2 - 2)) * Q_odd[j - K + 1L, j])
    }
  }
  
  b_even <- Q_even %*% a_even
  b_odd  <- Q_odd  %*% a_odd
  b <- rep(NA_real_, length(a))
  b[ SEQ] <- b_even
  b[-SEQ] <- b_odd

  # https://en.wikipedia.org/wiki/Binomial_theorem
  # b is the vector of coefficients for powers of (2 * p - 1)
  # b_[ , 1] is the vector of coefficients for powers of p
  # b_[ , 2] is the vector of coefficients for powers of p + 1
  # https://en.wikipedia.org/wiki/Budan%27s_theorem
  b_ <- matrix(0, length(b), 2)
  for (n in 0:(length(b) - 1L)) for (k in 0:n) {
    temp <- b[n + 1L] * choose(n, k) * 2^k
    b_[k + 1L, 1L] <- b_[k + 1L, 1L] + temp * (-1)^(n - k)
    b_[k + 1L, 2L] <- b_[k + 1L, 2L] + temp
  }
  if (diff(apply(b_, MARGIN = 2L, FUN = function(x) {
        sum(diff(ifelse(x == 0, NA_integer_, sign(x))) != 0, na.rm = TRUE)
      })) < 0 ) {
    roots <- polyroot(1:(length(b) - 1L) * b_[-1, 1L]) # local minima / maxima
    real_roots <- Re(roots[abs(Im(roots)) < 1e-13])
    real_roots <- real_roots[real_roots > 0 & real_roots < 1]
    if (length(real_roots) == 0L) return(TRUE)
    q <- Vectorize(function(p) (T(p, n = 0:(length(a) - 1L)) %*% a)[1])
    curve(q(p), from = 0, to = 1,
          n = 10001, xname = "p", ylab = expression(theta), axes = FALSE, las = 1)
    axis(1)
    abline(v = real_roots, col = 2, lty = 2)
    stop("\nImplied quantile function is decreasing near ", 
         paste(round(real_roots, digits = 4), collapse = ", "),  ".",
         "\nTry increasing the number of quantiles and / or changing their values.")
  }
  return(TRUE)
}
