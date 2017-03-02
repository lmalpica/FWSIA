---
output: word_document
---

We can write the model as follows:

$$
y_i \sim \mathrm{Normal} \left(y^\mathrm{true}_i, \tau_i^2\right),
$$

$$
y_i^\mathrm{true} = \beta_0 + \beta_1 H_i + \beta_2 L_i + \beta_3 R_i +
  \beta_4 R_i L_i + \beta_5 R_i H_i +
  \epsilon_{i},
$$

$$
\epsilon_i \sim \mathrm{Normal}(0, \sigma^2),
$$

where $y_i$ and $\tau_i$ represent the mean and standard deviation of the
log-transformed posterior of one of our response variables for site $i$,
$y^\mathrm{true}_i$ represents the "true" unobserved value of that response
variable, the $\beta$'s represent estimated coefficients, and $H_i$, $L_i$, and
$R_i$ represent the habitat index variable, the lionfish density, and a binary
variable representing whether lionfish were removed, respectively, for each site $i$.
The variable
$\epsilon_i$ represents independent normally distributed residual error with mean 0 and
standard deviation $\sigma$.

In order to make the coefficients approximately comparable across predictors, we centered and scaled each continuous predictor by subtracting its mean and dividing by two times its standard deviation (Gelman REF). We scaled by *two* times the standard deviation so that the coefficients on the continuous predictors are approximately comparable to the removal treatment binary predictor (Gelman REF). In the case of the binary removal treatment, we centered the predictor by subtracting its mean. Therefore, instead of $R_i$ being represented as 0s and 1s, it was represented as -0.5 (no removal) or 0.5 (lionfish removed). Because the removal treatment is centered, the coefficients for lionfish density and the habitat variable represent the effects of these variables averaged across the lionfish removal treatments --- they represent effects for a theoretical site with half the lionfish removal treatment.

We fit our models in a Bayesian framework using
Stan 2.14.1 [REFs] and R 3.3.2 [REF] to incorporate
the measurement uncertainty, $\tau$.
We assigned weakly informative
priors on all parameters: Normal$(0, 2^2)$ priors on the slope parameters
$\beta_1$ through $\beta_5$,
a Normal$(0, 10^2)$ prior on the intercept $\beta_0$, and
a half-\(t(3, 0, 2)\) prior (i.e. degrees of
freedom of 3 and scale of 2) on $\sigma$.
We ran four chains and 3,000 iterations, discarding the
first 1,500 iterations of each chain as
warm up.
We checked for chain convergence visually with trace plots, ensured
that \(\hat{R}<1.05\) (the potential scale reduction factor), and that the
effective sample size was greater than 200 for all parameters
(Gelman et al. 2014). We used a non-centered parameterization
of our hierarchical model and increased the target acceptance rate
in Stan to 0.99 to increase sampling efficiency
and ensure unbiased estimates (Stan manual REF, Monnahan et al. 2016).
Code for our model is available at... [supplement?]
