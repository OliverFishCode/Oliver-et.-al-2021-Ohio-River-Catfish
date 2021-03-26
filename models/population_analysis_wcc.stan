data {
  // VBGF and length-weight regression
  int<lower = 0> n1;             // N for VBGF and LW data
  vector<lower = 0>[n1] age;     // Age of fish for VBGF
  vector<lower = 0>[n1] length;  // Length of fish
  vector[n1] log10length;        // Log10 length
  vector[n1] log10weight;        // Log10 mass
  real s_vb;                     // Prior on VBGF error
  real s_lw;                     // Prior on LW error
  
  // Catch-curve analysis
  int<lower = 0> n2;             // N obs for CC data
  int<lower = 0> ngear;          // Number of gears for CC
  vector[n2] age2;               // Age of catch for CC
  vector[n2] N;                  // Number caught per age CC
  int gear[n2];                  // Gear for CC
  real s_cc;                     // Prior on CC error
  real s_cc_2;                    // prior on weighted CC error
  int<lower = 0> n_ages;         // Number of new ages for preds
  row_vector[n_ages] ages;           // New ages for prediction
}
parameters {
  // VBGF
  real<lower = 0> sigma_vb;      // Error
  real linf;                     // Mean asymptotic length, log_e
  real k;                        // Brody growth coeff, log_e
  real lt0;                      // Length at age zero, log_e
  
  // Length-weight regression
  real<lower = 0> sigma_lw;      // Error
  real a;                        // Intercept
  real b;                        // Slope
  
  // Catch-curve
  real<lower = 0> sigma_cc;      // Error
  vector[ngear] alpha;           // Intercept
  vector<upper = 0>[ngear] beta; // Slope
  
  real<lower = 0> sigma_cc_2;      // Weighted error
  vector[ngear] alpha_2;           // Weighted intercept
  vector<upper = 0>[ngear] beta_2; // Weighted slope 
}
transformed parameters{
  // VBGF
  real Linf = exp(linf);
  real K = exp(k);
  real t0 = exp(lt0) - 10;
  matrix[ngear, n_ages] N_2;
  matrix[ngear, n_ages] N_2_R;
  matrix[ngear, n_ages] sigma_weighted;
  
  // Weighted catch-curve prediction
  for(i in 1:n_ages){
   N_2[gear, i] = alpha[gear] + beta[gear] * ages[i];
  }
  
  N_2_R = exp(N_2);
  
  // Weighted variance
  sigma_weighted = sigma_cc_2 * N_2_R;
}
model {
  vector[n1] y;
  vector[n1] ew;
  vector[n2] c;
  vector[n2] c_2;
  
  // VBGF
  y = Linf * (1 - exp( -K * (age - t0))); // Expectation
  length ~ normal(y, sigma_vb);           // Likelihood
  sigma_vb ~ exponential(1/s_vb);
  linf ~ normal(0, 1);
  k ~ normal(0, 1);
  lt0 ~ normal(0, 1);
  
  // Length-weight regression
  ew = a + b*log10length;                // Expectation
  log10weight ~ normal(ew, sigma_lw);    // Likelihood
  sigma_lw ~ exponential(1/s_lw);
  a ~ normal(0, 1);
  b ~ normal(0, 1);
  
  // Catch-curve round 1
  c = alpha[gear] + beta[gear] .* age2;  // Expectation
  N ~ normal(c, sigma_cc);               // Likelihood
  sigma_cc ~ exponential(1/s_cc);
  alpha ~ normal(0, 1);
  beta ~ normal(0, 1);
  
  // Catch-curve round 2
  c_2 = alpha_2[gear] + beta_2[gear] .* age2;  // Expectation
  for(i in 1:n_ages){
   N ~ normal(c_2, sigma_weighted[gear, i]);             // Likelihood
  }
  sigma_cc_2 ~ exponential(1/s_cc_2);
  alpha_2 ~ normal(0, 1);
  beta_2 ~ normal(0, 1);
  
}
