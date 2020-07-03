functions {

  real no_name1_db_icdf(real p, vector c) { // Clenshaw?
    int Kp1 = rows(c);
    vector[Kp1] cheb;
    cheb[1] = 1;
    cheb[2] = 2 * p - 1;
    for (k in 3:Kp1) cheb[k] = (4 * p - 2) * cheb[k - 1] - cheb[k - 2];
    return dot_product(cheb, c);
  }
  
  real no_name1_lb_icdf(real p, vector c) {
    return exp(atanh(no_name1_db_icdf(p, c)));
  }
  
  real no_name1_ub_icdf(real p, vector c) {
    return atanh(no_name1_db_icdf(p, c));
  }
  
}
