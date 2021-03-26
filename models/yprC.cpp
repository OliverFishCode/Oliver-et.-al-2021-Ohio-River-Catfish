#include <Rcpp.h>
using namespace Rcpp;

// Code for running YPR models from Dippold et al. (2016)

// Main function declaration
//[[Rcpp::export]]
Rcpp::List yprC(NumericVector fm, NumericVector nm, NumericMatrix N,
                NumericMatrix Wt){

  // Define iterators for HMC samples (n) and age (f)
  int n = nm.size();
  int f = fm.size();

  // Declare objects created in fxn
  Rcpp::NumericMatrix Nd(n, f);
  Rcpp::NumericMatrix Nc(n, f);
  Rcpp::NumericVector Z(n);
  Rcpp::NumericMatrix YPR(n, f);

  for(int i = 0; i < n; i++){
    for(int t = 1; t < f; t++){
      // Total instantaneous mortality
      Z[i] = fm[t] + nm[i];
      // Population change
      N(i, t) = N(i, (t-1)) * exp(-Z[i]);
      // Number of individuals that died
      Nd(i, t) = N(i, t) * (1 - exp(-Z[i]));
      // Number that died from fishing
      Nc(i, t) = Nd(i, t) * (fm[t]/Z[i]);
      // YPR
      YPR(i, t) = Nc(i,t)*Wt(i, t);
    }
  }

  // Output list declaration
  List out;

  // Add output to list
  out["Z"] = Z;
  out["N"] = N;
  out["Nd"] = Nd;
  out["Nc"] = Nc;
  out["ypr"] = YPR;
  return out;
}
