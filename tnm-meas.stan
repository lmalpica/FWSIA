data {
  int N;               // number of observations
  int K;               // number of predictors
  real y_meas[N];      // measurement of y
  real<lower=0> tau[N];   // measurement sd on y
  matrix[N, K] X;      // model predictor matrix
  int N_new;               // number (new) of observations
  matrix[N_new, K] X_new;      // (new) model predictor matrix
}
parameters {
  vector[K] beta;      // vector of predictors
  real alpha;          // intercept
  real<lower=0> sigma; // residual sd
  real y_raw[N];
}
transformed parameters {
  real y[N];           // unknown true y value
  for (i in 1:N) {
    y[i] = alpha + X[i, ] * beta + sigma * y_raw[i]; // non-centered parameterization 'trick'
  }
} 
model { 
  sigma ~ student_t(3, 0, 2);  // prior
  alpha ~ normal(0, 10);  // prior
  beta ~ normal(0, 2);   // prior
  y_meas ~ normal(y, tau);     // measurement model
  y_raw ~ normal(0, 1);        // non-centered parameterization 'trick'
}
generated quantities{
  real y_pred[N];
  real y_pred_new[N_new];
  for (i in 1:N) {
    y_pred[i] = alpha + X[i, ] * beta;
  }
  for (i in 1:N_new) {
    y_pred_new[i] = alpha + X_new[i, ] * beta;
  }
}
