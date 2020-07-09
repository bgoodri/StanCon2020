functions {

  matrix make_T(data vector u) { // this is overpromoted
    int Kp1 = rows(u);
    vector[Kp1] constant = 4 * u - 2;
    matrix[Kp1, Kp1] T;
    if (Kp1 < 2) reject("u must have at least one element")
    T[ , 1] = rep_vector(1, Kp1);
    T[ , 2] = 0.5 * constant;
    if (u[1] < 0) reject("u must be non-negative");
    if (u[2] <= u[1]) reject("u must be increasing");
    for (k in 3:Kp1) {
      if (u[k] <= u[k - 1]) reject("u must be increasing");
      T[ , k] = constant .* T[ , k - 1] - T[ , k - 2];
    }
    if (u[Kp1] > 1) reject("all elements of u must be <= 1");
    return T;
  }

  real no_name1_db_icdf(real p, data vector u, data vector theta) {
    int Kp1 = rows(u);
    real constant = 4 * p - 2;
    vector[Kp1] cheb;
    cheb[1] = 1;
    cheb[2] = 0.5 * constant;
    for (k in 3:Kp1) cheb[k] = constant * cheb[k - 1] - cheb[k - 2];
    return dot_product(cheb, make_T(u[1:rows(theta)]) \ theta);
  }
  
  real no_name1_lb_icdf(real p, data vector u, data vector theta) {
    return exp(atanh(no_name1_db_icdf(p, u, tanh(log(theta)))));
  }
  
  real no_name1_ub_icdf(real p, data vector u, data vector theta) {
    return atanh(no_name1_db_icdf(p, u, tanh(theta)));
  }
  
}
