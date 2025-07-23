#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(Y);         
  DATA_VECTOR(Y_den);     
  DATA_MATRIX(D);                          // mapping of sample pairs
  DATA_MATRIX(X);                          // pair-level variables (diss. gradients)
  DATA_MATRIX(W);                          // sample-level variables (effects on uniqueness)
  DATA_VECTOR(weights)
  
  DATA_SPARSE_MATRIX(Z);                   // Random effects design
  DATA_INTEGER(has_random);                // Flag for re
  DATA_INTEGER(mono);                      // Flag (monotonic diss. effects)
  DATA_INTEGER(link);                      // link
  DATA_INTEGER(family);                    // family
  DATA_IVECTOR(map_re);                    // Mapping RE columns
  
  // Parameters
  PARAMETER(intercept);                   // Intercept
  PARAMETER_VECTOR(beta);                 // Slopes diss. grad
  PARAMETER_VECTOR(lambda);               // Slopes uniq.
  PARAMETER_VECTOR(u);                    // random effects
  PARAMETER_VECTOR(log_sigma_re);         // sigma of random effects
  PARAMETER(log_scale);                   // scale parameter

  // Derived quantities
  int n_pairs = Y.size();
  int n_X = X.cols();
  int n_W = W.cols();
  
  vector<Type> e_beta = beta;
  // Exp beta
  if (mono) {
    e_beta = exp(e_beta);
  }
  
  // Contributions 
  vector<Type> pair_comp_contrib(n_pairs);
  vector<Type> site_comp_contrib(n_pairs);
  vector<Type> re_comp(n_pairs);

  pair_comp_contrib.setZero();
  site_comp_contrib.setZero();
  re_comp.setZero();
  
  // Initiate log likelihood
  Type nll = 0.0;
  
  // Handle random effects
  vector<Type> sigma_re = exp(log_sigma_re);
  if (has_random) {
    
    // re nll
    for (int i = 0; i < u.size(); i++) {
      nll -= dnorm(u(i), Type(0), sigma_re(map_re[i]), true);
    }
    
    vector<Type> Z_cont = Z * u;
    
    // re full contribution
    for (int i = 0; i < n_pairs; i++) {
      int s1 = CppAD::Integer(D(i, 0));  
      int s2 = CppAD::Integer(D(i, 1));
      re_comp(i) = Z_cont(s1) + Z_cont(s2);
    }
  }
  
  // fixed effects contribution
  if (n_X > 0) {
    matrix<Type> pair_comp(n_pairs, n_X);
    for (int i = 0; i < n_pairs; i++) {
      int s1 = CppAD::Integer(D(i, 0));  
      int s2 = CppAD::Integer(D(i, 1));  
      pair_comp.row(i) = (X.row(s1) - X.row(s2)).array().abs();
    }
    pair_comp_contrib = pair_comp * e_beta;
  }
  
  if (n_W > 0) {
    matrix<Type> site_comp(n_pairs, n_W);
    for (int i = 0; i < n_pairs; i++) {
      int s1 = CppAD::Integer(D(i, 0));  
      int s2 = CppAD::Integer(D(i, 1));  
      site_comp.row(i) = (W.row(s1) + W.row(s2)).array();
    }
    site_comp_contrib = site_comp * lambda;
  }
  
  // linear predictor
  vector<Type> eta = intercept + pair_comp_contrib + site_comp_contrib + re_comp;
  
  
  // Expected dissimilarity
  vector<Type> mu(n_pairs); 
  
  // LINK 
  if (link == 0) {  // IDENTITY
    mu = eta;
  }
  
  if (link == 1) {  // LOGIT
    mu = exp(eta)/(1+exp(eta));
  }
  
  if (link == 2) {  // GDM
    mu = Type(1) - exp(-eta);
  }
  
  if (link == 3) { // GAUSS
    mu = Type(1) - exp(-(eta*eta));
  }

  // family nll
  if (family == 0) {  // Gaussian
    Type scale = exp(log_scale);
    for (int i = 0; i < n_pairs; i++) {
      nll -= weights(i)*dnorm(Y(i), mu(i), scale, true);
    }
  }
  
  else if (family == 1) {  // Binomial
    for (int i = 0; i < n_pairs; i++) {
     nll -= weights(i)*dbinom(Y(i), Y_den(i), mu(i), true);
    }
  }
  
  else if (family == 2) {  // Beta
    Type scale = exp(log_scale);
    // reparam. alpha beta into mu phi (scale)
    vector<Type> a = mu * scale;
    vector<Type> b = (Type(1) - mu) * scale;
    // compute likelihood
    for (int i = 0; i < n_pairs; i++) {
      nll -= weights(i)*dbeta(Y(i), a(i), b(i), true);
    }
  }
  
  ADREPORT(e_beta);
  ADREPORT(lambda);
  ADREPORT(intercept);
  
  REPORT(u);
  REPORT(e_beta);
  REPORT(lambda);
  REPORT(intercept);
  REPORT(sigma_re);
  return nll;
}



