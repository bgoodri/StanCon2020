# assumes you have already done rstan::expose_stan_functions("quantile_functions.stan")
qgld <- Vectorize(GLD_icdf, vectorize.args = "p")

GLD_solver <- function(lower_quartile, median, upper_quartile, 
                       other_quantile, alpha, check = TRUE) {
  obj_xi <- function(xi) {   # chi will be local
    GLD_icdf(alpha, median, IQR = upper_quartile - lower_quartile, chi, xi) - other_quantile
  }
  obj_chi <- function(chi) { # xi will be local (as will skewness)
    S_75 <- S_(0.75, chi, xi)
    S_25 <- S_(0.25, chi, xi)
    (S_75 + S_25 - 2 * S_(0.5, chi, xi) ) / (S_75 - S_25) - skewness
  }
  skewness <- (upper_quartile + lower_quartile - 2 * median) / 
              (upper_quartile - lower_quartile)
  chi <- skewness # good starting value if skewness is close to zero
  names(chi) <- NULL
  old_chi <- 2
  while (abs(chi - old_chi) > 1e-8) {
    xi <- try(uniroot(obj_xi, lower = 1e-2, upper = 1 - 1e-2)$root, silent = TRUE)
    if (!is.numeric(xi)) {
      if (check) stop("no GLD is possible with these quantiles")
      return(c(asymmetry = NA_real_, steepness = NA_real_))
    }
    old_chi <- chi
    chi <- try(uniroot(obj_chi, lower = -1 + 1e-6, upper = 1 - 1e-6)$root, silent = TRUE)
    if (!is.numeric(chi)) {
      if (check) stop("no GLD is possible with these quantiles")
      return(c(asymmetry = NA_real_, steepness = NA_real_))
    }
  }
  if (alpha == 0 || alpha == 1) xi <- xi - .Machine$double.eps
  a_s <- c(asymmetry = chi, steepness = xi)
  if (check) {
    low <- GLD_icdf(0, median, IQR = upper_quartile - lower_quartile,
                  asymmetry = a_s[1], steepness = a_s[2])
    if (alpha > 0 && low > -Inf) warning("solution implies a bounded lower tail at ", low)
    high <- GLD_icdf(1, median, IQR = upper_quartile - lower_quartile, 
                     asymmetry = a_s[1], steepness = a_s[2])
    if (alpha < 1 && high < Inf) warning("solution implies a bounded upper tail at ", high)
  }
  return(a_s)
}

GLD_solver_LBFGS <- function(lower_quartile, median, upper_quartile,
                             other_quantile, alpha, check = TRUE) {
  if (alpha == 0) {
    theta <- c(asymmetry =  0.5, steepness = 0.25)
  } else if (alpha == 1) {
    theta <- c(asymmetry = -0.5, steepness = 0.25)
  } else theta <- c(asymmetry = 0, steepness = 0.5 - 1e-8)
  IQR <- upper_quartile - lower_quartile
  skewness <- (upper_quartile + lower_quartile - 2 * median) / IQR
  fn <- function(theta) {
    S_75 <- S_(0.75, theta[1], theta[2])
    S_25 <- S_(0.25, theta[1], theta[2])
    r <- c((GLD_icdf(alpha, median, IQR, theta[1], theta[2]) - other_quantile) / IQR, 
           (S_75 + S_25 - 2 * S_(0.5, theta[1], theta[2]) ) / (S_75 - S_25) - skewness)
    if (any(!is.finite(r))) return(.Machine$double.xmax)
    return(crossprod(r)[1])
  }
  buffer <- .001
  opt <- optim(theta, fn = fn, method = "L-BFGS-B", 
               lower = c(-1, 0) + buffer, upper = c(1, 1) - buffer, 
               control = list(trace = 10 * check))
  if (check && opt$convergence != 0) 
    warning(opt$message,
            "\nno GLD is exactly possible with these quantiles, so check that the results are decent")
  return(opt$par)
}
  
GLD_solver_bounded <- function(bounds, median, IQR, check = TRUE, ...) {
  obj <- function(theta) {
    chi <- theta[1]
    xi  <- theta[2]
    if (abs(chi) > 1 || xi < 0 || xi > 1 || 
        xi > 0.5 * (1 + chi) || xi > 0.5 * (1 - chi)) return(NA_real_)
    out <- (bounds[1] - GLD_icdf(0, median, IQR, chi, xi)) ^ 2 +
           (bounds[2] - GLD_icdf(1, median, IQR, chi, xi)) ^ 2
    return(out)
  }
  start <- c(0, 0.5 - 1 / sqrt(5))
  # start is a solution that implies the GLD is uniform(a, b) where 
  # median == 0.5 * (a + b) and IQR == 0.5 * (b - a)
  sln <- suppressWarnings(nlm(obj, p = start, fscale = 0, ...))
  if (check && sln$code >= 4) {
    warning("solution for asymmetry and steepness is dubious",
            " try passing optional arguments to nlm via ... ")
  }
  a_s <- c(asymmetry = sln$estimate[1], steepness = sln$estimate[2])
  if (check && sln$minimum > 1e-15) {
    low  <- GLD_icdf(0, median, IQR, asymmetry = a_s[1], steepness = a_s[2])
    high <- GLD_icdf(1, median, IQR, asymmetry = a_s[1], steepness = a_s[2])
    warning("no asymmetry and steepness values achieve the bounds exactly\n",
            paste("effective bounds are", low, "and", high))
  }
  return(a_s)
}

GLD_solver_maxent <- function(lower_quartile, median, upper_quartile, 
                              check = TRUE, ...) {
  IQR <- upper_quartile - lower_quartile
  skewness <- (upper_quartile + lower_quartile - 2 * median) / IQR
  obj_chi <- function(chi, xi) {
    S_75 <- S_(0.75, chi, xi)
    S_25 <- S_(0.25, chi, xi)
    (S_75 + S_25 - 2 * S_(0.5, chi, xi) ) / (S_75 - S_25) - skewness
  }
  obj_xi <- function(xi) {
    chi <- try(uniroot(obj_chi, lower = -1 + 1e-6, upper = 1 - 1e-6, xi = xi)$root, 
               silent = TRUE)
    if (!is.numeric(chi)) stop("no GLD is possible with these quantiles")
    GLD_de(IQR, asymmetry = chi, steepness = xi)
  }
  opt <- optimize(obj_xi, interval = 0:1, maximum = TRUE)
  xi <- opt$maximum
  chi <- uniroot(obj_chi, lower = -1 + 1e-6, upper = 1 - 1e-6, xi = xi)$root
  return(c(asymmetry = chi, steepness = xi))
}

GLD_solver_mode <- function(lower_quartile, median, upper_quartile, 
                            mode, check = TRUE, ...) {
  IQR <- upper_quartile - lower_quartile
  skewness <- (upper_quartile + lower_quartile - 2 * median) / IQR
  obj_chi <- function(chi, xi) {
    S_75 <- S_(0.75, chi, xi)
    S_25 <- S_(0.25, chi, xi)
    (S_75 + S_25 - 2 * S_(0.5, chi, xi) ) / (S_75 - S_25) - skewness
  }
  obj_xi <- function(xi) {
    chi <- try(uniroot(obj_chi, lower = -1 + 1e-6, upper = 1 - 1e-6, xi = xi)$root, 
               silent = TRUE)
    if (!is.numeric(chi)) stop("no GLD is possible with these quantiles and mode")
    alpha <- 0.5 * (0.5 - xi) / sqrt(xi * (1 - xi))
    beta  <- 0.5 * chi / sqrt(1 - chi ^ 2)
    lambda_3 <- alpha + beta
    lambda_4 <- alpha - beta
    if ( ( (lambda_3 > 2) && (lambda_4 > 2) ) ||
         ( (lambda_3 < 1) && (lambda_4 < 1) ) ) {
      lambda_2 <- (S_(0.75, chi, xi) - S_(0.25, chi, xi)) / IQR
      lambda_1 <- median - S_(0.5, chi, xi) / lambda_2
      fBasics::gldMode(lambda_1, lambda_2, lambda_3, lambda_4) - mode
    } else stop("no GLD is possible with these quantiles and mode")
  }
  sln <- uniroot(obj_xi, lower = 0, upper = (17 - 4 * sqrt(17)) / 34)
  
  xi <- opt$maximum
  chi <- uniroot(obj_chi, lower = -1 + 1e-6, upper = 1 - 1e-6, xi = xi)$root
  return(c(asymmetry = chi, steepness = xi))
}
