functions {
  /*
     Checks whether real numbers are finite and ordered
     
     @param theta real array of numbers
     @throws if real numbers are not finite or not ordered
     @return 1 if real numbers are finite and ordered
   */
  int in_order(real[] theta) {
    if (num_elements(theta) != 5) reject("wrong number of elements");
    if (theta[1] == negative_infinity()) reject("first element must be finite");
    for (k in 2:5) if (theta[k] <= theta[k - 1]) 
      reject("bounds and quantiles are not in the right order");
    return 1;
  }

  /* Johnson Quantile Parameterized Distributions (J-QPD) */

  /*
     Inverse CDF of the J-QPD-S semi-bounded distribution, which has moments
     
     See equation 9 of
     http://metalogdistributions.com/images/Johnson_Quantile-Parameterized_Distributions.pdf
     
     @param p real cumulative probability
     @param alpha fixed proportion of distribution below quantiles[1]
     @param lower_bound real lower bound to the random variable
     @param quantiles vector of size three ordered quantiles
     @return real number greater than lower_bound
   */
  real JQPDS_icdf(real p, real lower_bound, data real alpha, vector quantiles) {
    if (p < 0 || p > 1)         reject("p must be between 0 and 1");
    if (alpha < 0 || alpha > 1) reject("alpha must be between 0 and 1");
    if (rows(quantiles) != 3)   reject("quantiles must have three elements");
    if (in_order({lower_bound, quantiles[1], quantiles[2], 
                  quantiles[3], positive_infinity()})) {
      real c = inv_Phi(1 - alpha);
      vector[3] quantiles_ = quantiles - lower_bound;
      real L = log(quantiles_[1]);
      real B = log(quantiles_[2]);
      real H = log(quantiles_[3]);
      real HmL = H - L;
      real denom = fmin(B - L, H - B);
      real numer = sinh(acosh(0.5 * HmL / denom));
      real delta = numer / c;
      real lambda = denom / numer;
      real LHm2B = L + H - 2 * B;
      real k = sqrt(1 + square(numer));
      real n;
      real theta;
      real z;
      if (LHm2B < 0) {
        n = -1;
        theta = quantiles_[3];
        z = inv_Phi(p);
      } else if (LHm2B > 0) {
        n = 1;
        theta = quantiles_[1];
        z = inv_Phi(p);
      } else { // LHm2B = 0 -> removable discontinuity
        real sigma = delta != 0 ? lambda * delta : (H - B) / c;
        theta = quantiles_[2];
        return lower_bound + theta * exp(sigma * inv_Phi(p));
      }
      // return lower_bound + theta * exp(lambda * sinh(asinh(delta * z) + asinh(n * numer)));
      // http://www.metalogdistributions.com/images/J-QPD_Parameterizations.pdf section 1.2
      return lower_bound + theta * 
        exp(lambda * delta * (k * z + n * c * sqrt(1 + square(delta * z))));
    }
    return not_a_number(); // never reaches
  }

  /*
     Pseudo-random number generator of the J-QPD-S semi-bounded distribution
     
     @param alpha fixed proportion of distribution below quantiles[1]
     @param lower_bound real lower bound to the random variable
     @param quantiles vector of size three ordered quantiles
     @return real number greater than lower_bound
   */
  real JQPDS_rng(real lower_bound, data real alpha, vector quantiles) {
    return JQPDS_icdf(uniform_rng(0, 1), lower_bound, alpha, quantiles);
  }

  /*
     Inverse CDF of the J-QPD-S-II semi-bounded distribution, which lacks moments
     
     See equation 14 of
     http://metalogdistributions.com/images/Johnson_Quantile-Parameterized_Distributions.pdf
     It is the limit of the J-QPD-B distribution as the upper bound diverges.
     
     @param p real cumulative probability
     @param alpha fixed proportion of distribution below quantiles[1]
     @param lower_bound real lower bound to the random variable
     @param quantiles vector of size three ordered quantiles
     @return real number greater than lower_bound
   */
  real JQPDS2_icdf(real p, real lower_bound, data real alpha, vector quantiles) {
    if (p < 0 || p > 1)         reject("p must be between 0 and 1");
    if (alpha < 0 || alpha > 1) reject("alpha must be between 0 and 1");
    if (rows(quantiles) != 3)   reject("quantiles must have three elements");
    if (in_order({lower_bound, quantiles[1], quantiles[2], 
                  quantiles[3], positive_infinity()})) {
      real c = inv_Phi(1 - alpha);
      vector[3] quantiles_ = quantiles - lower_bound;
      real L = log(quantiles_[1]);
      real B = log(quantiles_[2]);
      real H = log(quantiles_[3]);
      real HmL = H - L;
      real denom = fmin(B - L, H - B);
      real temp = acosh(0.5 * HmL / denom);
      real delta = temp * inv(c);
      real lambda = denom * inv(sinh(temp));
      real LHm2B = L + H - 2 * B;
      real n;
      real theta;
      if (LHm2B < 0) {
        n = -1;
        theta = quantiles_[3];
      } else if (LHm2B > 0) {
        n = 1;
        theta = quantiles_[1];
      } else { // LHm2B = 0 -> removable discontinuity
        return lower_bound + quantiles[2] * exp(lambda * sinh(delta * inv_Phi(p)));
      }
      return lower_bound + theta * exp(lambda * sinh(delta * (inv_Phi(p) + n * c)));
    }
    return not_a_number(); // never reaches
  }

  /*
     Pseudo-random number generator of the J-QPD-S-II semi-bounded distribution
     
     @param alpha fixed proportion of distribution below quantiles[1]
     @param lower_bound real lower bound to the random variable
     @param quantiles vector of size three ordered quantiles
     @return real number greater than lower_bound
   */
  real JQPDS2_rng(real lower_bound, data real alpha, vector quantiles) {
    return JQPDS2_icdf(uniform_rng(0, 1), lower_bound, alpha, quantiles);
  }

  /*
     Inverse CDF of the J-QPD-B bounded distribution
     
     See equation 7 of
     http://metalogdistributions.com/images/Johnson_Quantile-Parameterized_Distributions.pdf
     
     @param p real cumulative probability
     @param alpha fixed proportion of distribution below quantiles[1]
     @param bounds vector of size two containing the lower and upper bounds
     @param quantiles vector of size three ordered quantiles
     @return real number greater than lower_bound
   */
  real JQPDB_icdf(real p, row_vector bounds, data real alpha, vector quantiles) {
    if (cols(bounds) != 2)      reject("bounds must have two elements");
    if (bounds[2] == positive_infinity()) {
      return JQPDS2_icdf(p, alpha, bounds[1], quantiles);
    }
    if (p < 0 || p > 1)         reject("p must be between 0 and 1");
    if (alpha < 0 || alpha > 1) reject("alpha must be between 0 and 1");
    if (rows(quantiles) != 3)   reject("quantiles must have three elements");
    if (in_order({bounds[1], quantiles[1], quantiles[2], quantiles[3], bounds[2]})) {
      real c = inv_Phi(1 - alpha);
      real l = bounds[1];
      real u = bounds[2];
      real uml = u - l;
      real L = inv_Phi( (quantiles[1] - l) / uml );
      real B = inv_Phi( (quantiles[2] - l) / uml );
      real H = inv_Phi( (quantiles[3] - l) / uml );
      real HmL = H - L;
      real delta = acosh(0.5 * HmL / fmin(B - L, H - B)) / c;
      real lambda = HmL / sinh(2 * delta * c);
      real LHm2B = L + H - 2 * B;
      real n;
      real zeta;
      if (LHm2B < 0) {
        n = -1;
        zeta = H;
      } else if (LHm2B > 0) {
        n = 1;
        zeta = L;
      } else { // LHm2B = 0 -> removable discontinuity
        return l + uml * Phi(B + 0.5 * HmL / c * inv_Phi(p));
      }
      return l + uml * Phi(zeta + lambda * sinh(delta * (inv_Phi(p) + n * c)));
    }
    return not_a_number(); // never reached
  }

  /*
     Pseudo-random number generator of the J-QPD-B bounded distribution
     
     @param alpha fixed proportion of distribution below quantiles[1]
     @param bounds vector of size two containing the lower and upper bounds
     @param quantiles vector of size three ordered quantiles
     @return real number between the two elements of bounds
   */
  real JQPDB_rng(data row_vector bounds, real alpha, vector quantiles) {
    return JQPDB_icdf(uniform_rng(0, 1), bounds, alpha, quantiles);
  }
  
  /* Quantile Parameteried Normal (qnormal) distribution */
  
  /*
     Inverse CDF of the qnormal distribution. See
  
     http://metalogdistributions.com/images/KeelinPowley_QuantileParameterizedDistributions_2011.pdf
     
     @param p real cumulative probability
     @param a vector of size four coefficients
     @return real number
   */
  real qnormal_icdf(real p, vector a) {
    if (p < 0 || p > 1) reject("p must be between 0 and 1");
    if (rows(a) == 4) {
      real a3 = a[3];
      real a4 = a[4];
      real mu = a[1] + a4 * p;
      real sigma = a[2] + a3 * p;
      real z = inv_Phi(p);
      real denom = sigma + exp(std_normal_lpdf(z)) * (a3 * z + a4);
      if (denom < 0) reject("a does not imply a valid distribution");
      return mu + sigma * z;
    } else reject("a must be of size four");
    return not_a_number(); // never reached
  }
  
  /*
     Pseudo-random number generator of the qnormal distribution. See
  
     @param a vector of size four coefficients
     @return real number
   */
  real qnormal_rng(vector a) {
    return qnormal_icdf(uniform_rng(0, 1), a);
  }

  vector qnormal_coefficients(data vector quantiles, data vector p) {
    if (rows(quantiles) != 4) reject("quantiles must be of size four");
    if (rows(p) != 4)         reject("p must be of size four");
    if (in_order({quantiles[1], quantiles[2], quantiles[3], 
                  quantiles[4], positive_infinity()}) == 
        in_order({p[1], p[2], p[3], p[4], 1})) {
      vector[4] inv_CDF = inv_Phi(p);
      return [ [1, p[1], inv_CDF[1], p[1] * inv_CDF[1]], 
               [1, p[2], inv_CDF[2], p[2] * inv_CDF[2]],
               [1, p[3], inv_CDF[3], p[3] * inv_CDF[3]],
               [1, p[4], inv_CDF[4], p[4] * inv_CDF[4]] ] \ quantiles;
    }
    return quantiles; // never reaches
  }
  
  /* Metalog Distribution */
  
  /* Inverse CDF of the three term metalog distribution
  
     See equations 18, 19, and 20 of 
     http://www.metalogdistributions.com/images/TheMetalogDistributions.pdf
    
     @param p real cumulative probability
     @param alpha fixed proportion of distribution below quantiles[1]
     @param quantiles vector of size three ordered quantiles
   */  
  real metalog3_icdf(real p, data real alpha, data vector quantiles) {
    real out_diff = quantiles[3] - quantiles[1];
    real log_odds_alpha = logit(alpha); 
    real a2 = 0.5 / log_odds_alpha * out_diff;
    real in_diff = quantiles[2] - quantiles[1];
    real a3 = (out_diff - 2 * in_diff) / ((1 - 2 * alpha) * log_odds_alpha);
    real k = 0.5 * (1 - 1.66711 * (0.5 - alpha));
    real r = in_diff / out_diff;
    if (alpha >= 0.5) reject("alpha must be less than 0.5");
    if (in_diff <= 0) reject("quantiles[2] must be greater than quantiles[1]");
    if (quantiles[2] >= quantiles[3]) 
      reject("quantiles[3] must greater than quantiles[2]");
    if (k >= r || r >= (1 - k)) 
      reject("quantiles do not imply a valid distribution");
    return quantiles[2] + (a2 + a3 * (p - 0.5)) * logit(p);
  }
  
  vector metalog_coefficients(data vector p, data vector quantiles) {
    int n = rows(p);
    matrix[n, n] Y;
    if (n == 0) reject("p cannot be of size zero");
    if (rows(quantiles) == n) {
      vector[n] log_odds = n > 1 ? logit(p) : rep_vector(not_a_number(), n);
      vector[n] pmhalf =   n > 2 ? p - 0.5  : rep_vector(not_a_number(), n);
      int odd = 1;
      Y[ , 1] = rep_vector(1, n);
      if (n > 1) Y[ , 2] = log_odds;
      if (n > 2) Y[ , 3] = pmhalf .* log_odds;
      if (n > 3) Y[ , 4] = pmhalf;
      for (m in 5:n) {
        if (odd) {
          pmhalf .*= pmhalf;
          Y[ , m]  = pmhalf;
        } else Y[ , m] = pmhalf .* log_odds;
        odd = odd == 0;
      }
    } else reject("p and quantiles must be of the same size");
    return Y \ quantiles;
  }
  
  /* Generalized Lambda Distribution (GLD) with good parameterization
  
  /*
     Helper function for the inverse CDF of the GLD
  
     See equation 11 of
     https://mpra.ub.uni-muenchen.de/37814/1/MPRA_paper_37814.pdf
  
     @param p real cumulative probability
     @param chi real skewness parameter between -1 and 1
     @param xi real steepness parameter between  0 and 1
     @return real number that is scaled and shifted to ~GLD
   */
  real S_(real p, real chi, real xi) {
    real alpha = 0.5 * (0.5 - xi) * inv_sqrt(xi * (1 - xi));
    real beta  = 0.5 * chi * inv_sqrt(1 - square(chi));
    if (p > 0 && p < 1) {
      if (chi != 0 || xi != 0.5) {
        if (fabs(alpha) != beta) {
          real s = alpha + beta;
          real d = alpha - beta;
          if (alpha == negative_infinity()) return 0;
          return (p ^ s - 1) / s - ( (1 - p) ^ d - 1 ) / d;
        } else if (xi == 0.5 * (1 + chi)) {
          real d = 2 * alpha;
          return log(p) - ( (1 - p) ^ d - 1 ) / d;
        } else {// xi == 0.5 * (1 - chi)
          real s = 2 * alpha;
          return (p ^ s - 1) / s - log1m(p);
        }
      } else return log(p) - log1m(p); // chi == 0 and xi == 0.5
    } else if (p == 0) { // equation 13
      return xi < 0.5 * (1 + chi) ? -inv(alpha + beta) : negative_infinity();
    } else if (p == 1) { // equation 14
      return xi < 0.5 * (1 - chi) ?  inv(alpha - beta) : positive_infinity();
    } else reject("p must be between zero and one");
    return not_a_number(); // never reaches
  }

  /*
     Inverse CDF of the GLD
  
     See equation 12 of
     https://mpra.ub.uni-muenchen.de/37814/1/MPRA_paper_37814.pdf
  
     @param p real cumulative probability
     @param median real median of the GLD
     @param IQR real inter-quartile range greater than 0
     @param asymmetry real parameter between -1 and 1
     @param steepness real parameter between  0 and 1
     @return real number that is ~GLD
   */
  real GLD_icdf(real p, real median, real IQR, real asymmetry, real steepness) {
    real CHI = fabs(asymmetry);
    if (IQR < 0) reject("IQR must be non-negative");
    if (steepness < 0 || steepness > 1) 
      reject("steepness must be between 0 and 1");
    if (CHI < 1)
      return median + IQR * 
        (S_(p,    asymmetry, steepness) - S_(0.50, asymmetry, steepness)) / 
        (S_(0.75, asymmetry, steepness) - S_(0.25, asymmetry, steepness));
    if (CHI > 1) reject("asymmetry must be between -1 and 1");
    if (steepness != 0)
      reject("steepness must be 0 when asymmetry is ", asymmetry);
    if (asymmetry == -1) return median + IQR * (  log(p) + log(2.0)) / log(3.0);
    if (asymmetry ==  1) return median - IQR * (log1m(p) + log(2.0)) / log(3.0);
    // if asymmetry == 1 GLD is an exponential with mean = median / log(2)
    return not_a_number(); // never reaches
  }

  /*
     Pseudo-random number generator of the GLD
  
     @param median real median of the GLD
     @param IQR real inter-quartile range greater than 0
     @param asymmetry real parameter between -1 and 1
     @param steepness real parameter between  0 and 1
     @return real number that is ~GLD
   */
  real GLD_rng(real median, real IQR, real asymmetry, real steepness) {
    return GLD_icdf(uniform_rng(0, 1), median, IQR, asymmetry, steepness);
  }

  /*
     Log quantile density function of the GLD (which does not depend on the median)
  
     See equation 14 of
     https://mpra.ub.uni-muenchen.de/37814/1/MPRA_paper_37814.pdf
  
     @param p real cumulative probability
     @param p_ real complementary probability
     @param theta real array containing the IQR, asymmetry, and steepness
     @param x_r real array that is not used
     @param x_i integer array that is not used
     @return real number
   */
  real GLD_qd_integrand(real p, real p_, real[] theta, real[] x_r, int[] x_i) {
    real asymmetry = theta[2];
    real steepness = theta[3];
    real alpha = 0.5 * (0.5 - steepness) * inv_sqrt(steepness * (1 - steepness));
    real beta  = 0.5 * asymmetry * inv_sqrt(1 - square(asymmetry));
    real log_q = log(theta[1]) 
               - log(S_(0.75, asymmetry, steepness) - S_(0.25, asymmetry, steepness));
    if (theta[1] < 0) reject("IQR must be non-negative");
    if (fabs(asymmetry) > 1) reject("asymmetry must be between -1 and 1");
    if (steepness < 0 || steepness > 1) reject("steepness must be between 0 and 1");
    if (asymmetry != 0 || steepness != 0.5) {
      if (fabs(alpha) != beta) {
        real s = alpha + beta;
        real d = alpha - beta;
        log_q += log_sum_exp((s - 1) * log(p),
                             (d - 1) * (p <= 0.5 ? log1m(p) : log(p_)));
      } else if (steepness == 0.5 * (1 + asymmetry)) {
        log_q += log_sum_exp(-log(p), 
                             (2 * alpha - 1) * (p <= 0.5 ? log1m(p) : log(p_)));
      } else {// xi == 0.5 * (1 - chi)
        log_q += log_sum_exp((2 * beta - 1) * log(p),
                             p <= 0.5 ? -log1m(p) : -log(p_));
      }
    } else log_q += log_sum_exp(-log(p), (p <= 0.5 ? -log1m(p) : -log(p_)));
    return log_q;
  }

  /*
     Differential entropy of the GLD (which does not depend on the median)
  
     See equation 14 of
     https://mpra.ub.uni-muenchen.de/37814/1/MPRA_paper_37814.pdf
  
     @param IQR real inter-quartile range greater than 0
     @param asymmetry real parameter between -1 and 1
     @param steepness real parameter between  0 and 1
     @return real number
   */
  real GLD_de(real IQR, real asymmetry, real steepness) {
    int x_i[0];
    return integrate_1d(GLD_qd_integrand, 0, 1, {IQR, asymmetry, steepness}, 
                        rep_array(0.0, 0), x_i, 1e-8);
  }

  /*
     System of equations to find shape parameters of the GLD
  
     See equations 16a and 16b of
     https://mpra.ub.uni-muenchen.de/43333/3/MPRA_paper_43333.pdf
  
     @param free vector of size two pre-parameters of chi and xi
     @param theta vector of size zero that is unused
     @param x_r real array of size seven containing quantiles at the
       follwing alpha levels: {1 / 8, 2 / 8, 3 / 8, 4 / 8, 5 / 8, 6 / 8, 7 / 8}
     @param x_i integer array of size one with a flag indicating whether to
       force the solution to imply unbounded support
     @return vector of two residuals
   */
  vector equations(vector free, vector theta, data real[] x_r, data int[] x_i) {
    real chi = tanh(free[1]);
    real lower_bound = x_i[1] ? 0.5 * (1 - fabs(chi)) : 0.0;
    real xi = x_i[1] ? lower_bound + inv_logit(free[2]) * (1 - lower_bound)
                     : inv_logit(free[2]);
    real high = S_(0.75, chi, xi);
    real low  = S_(0.25, chi, xi);
    real denom = high - low;
    // avoid overpromotion but conceptually ...
    // real IQR  = x_r[6] - x_r[2];
    // real skewness = (x_r[6] + x_r[2] - 2 * x_r[4]) / IQR;
    // real kurtosis = (x_r[7] - x_r[5] + x_r[3] - x_r[1]) / IQR;
    return [(x_r[6] + x_r[2] - 2 * x_r[4]) / (x_r[6] - x_r[2])
            - (high + low - 2 * S_(0.5, chi, xi)) / denom,
            (x_r[7] - x_r[5] + x_r[3] - x_r[1]) / (x_r[6] - x_r[2])
            - (S_(7.0 / 8, chi, xi) - S_(5.0 / 8, chi, xi) + 
               S_(3.0 / 8, chi, xi) - S_(1.0 / 8, chi, xi)) / denom]';
  }

  /*
     System of equations to find shape parameters of the GLD
  
     One quantile and equation 16a of
     https://mpra.ub.uni-muenchen.de/43333/3/MPRA_paper_43333.pdf
  
     @param free vector of size two pre-parameters of chi and xi
     @param theta vector of size zero that is unused
     @param x_r real array of size five containing:
       {lower_quartile, median, upper_quartile, other_quantile, alpha}
     @param x_i integer array of size one with a flag indicating whether to
       force the solution to imply unbounded support
     @return vector of two residuals
   */
  vector equations2(vector free, vector theta, data real[] x_r, data int[] x_i) {
    real chi = tanh(free[1]);
    real lower_bound = x_i[1] ? 0.5 * (1 + fabs(chi)) : 0.0;
    real xi = x_i[1] ? lower_bound + inv_logit(free[2]) * (1 - lower_bound)
                     : inv_logit(free[2]);
    real high = S_(0.75, chi, xi);
    real low  = S_(0.25, chi, xi);
    real denom = high - low;
    return [(x_r[3] + x_r[1] - 2 * x_r[2]) / (x_r[3] - x_r[1]) 
             - (high + low - 2 * S_(0.5, chi, xi)) / denom,
            x_r[4] - GLD_icdf(x_r[5], x_r[2], x_r[3] - x_r[1], chi, xi)]';
  }

  /*
     Solve system of equations to find shape parameters of the GLD
  
     @param x_r real array of size 5 or 7
     @param unbounded int indicating whether to enforce unbounded support
     @return vector of shape parameters chi and xi of the GLD
   */
  vector find_chi_xi(data real[] x_r, int unbounded) {
    vector[2] free;
    real chi;
    real lower_bound;
    real xi;
    if (num_elements(x_r) == 5) {
      free = algebra_solver(equations2, 
                            [sqrt(machine_precision()), -machine_precision()]',
                            rep_vector(0, 0), x_r, {unbounded});
                            
    } else if (num_elements(x_r) == 7) {
      free = algebra_solver(equations,
                            [sqrt(machine_precision()), -machine_precision()]',
                            rep_vector(0, 0), x_r, {unbounded});
    } else reject("x_r must be of size 5 or 7");
    chi = tanh(free[1]);
    lower_bound = unbounded ? 0.5 * (1 + fabs(chi)) : 0.0;
    xi = lower_bound + inv_logit(free[2]) * (1 - lower_bound);
    return [chi, xi]';
  }
  
  real count_moments(real chi, real xi) {
    real alpha = 0.5 * (0.5 - xi) * inv_sqrt(xi * (1 - xi));
    real beta  = 0.5 * chi * inv_sqrt(1 - square(chi));
    real threshold = fmin(alpha + beta, alpha - beta);
    if (threshold >= 0) return positive_infinity();
    return -inv(threshold);
  }
  
  // student t https://pdfs.semanticscholar.org/5953/0b61f8e14530169dfc3a9149488f90be8fe0.pdf
  real student_t_icdf(real p, real nu, real mu, real sigma) {
    real nup1 = nu + 1;
    real nu2 = square(nu);
    real nu3 = nu2 * nu;
    real nu4 = square(nu2);
    real nu5 = nu4 * nu;
    real base = exp(0.5 * log(nu * pi()) + lgamma(0.5 * nu) - lgamma(0.5 * nup1)) * (p - 0.5);
    real out = base 
             + nup1 / (6 * nu) * base ^ 3 
             + nup1 * (7 * nu + 1) / (120 * nu2) * base ^ 5
             + nup1 * (127 * nu2 + 8 * nu + 1) / (5040 * nu3) * base ^ 7
             + nup1 * (4369 * nu3 - 537 * nu2 + 135 * nu + 1) / (362880 * nu4) * base ^ 9
             + nup1 * (243649 * nu4 - 90488 * nu3 + 26238 * nu2 - 2504 * nu + 1) / 
               (39916800 * nu5) * base ^ 11;
    return mu + sigma * out;
  }
}

