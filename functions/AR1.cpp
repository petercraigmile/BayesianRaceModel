
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double AR1LogLikelihood (NumericVector z,
			 double phi,
			 double sigma) {

  
  int n = z.size();
  double ll;
  
  ll = R::dnorm(z[0], 0.0, sigma/sqrt(1.0-phi*phi), TRUE);
  
  for (int i = 1; i < n; i++) {
    ll += R::dnorm(z[i], phi*z[i-1], sigma, TRUE);
  }
  
  return ll;
}

// [[Rcpp::export]]
NumericVector AR1sim (int n, double phi, double sigma) {

  NumericVector z(n);

  z[0] = R::rnorm(0.0, sigma/sqrt(1.0-phi*phi));
  
  for (int i = 1; i < n; i++) {
    z[i] = R::rnorm(phi*z[i-1], sigma);
  }

  return z;
}
