TNM Stan models
================
2016-12-10

``` r
library(dplyr)
library(ggplot2)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
ctl <- list(adapt_delta = 0.99) # more robust than the default 
library(bayesplot)
# devtools::install_github("seananderson/stanhelpers") # if needed
library(stanhelpers)
```

Data
====

This loads a file containing posterior distributions for six Trophic Niche Metrics (TNM): dNr, dCr, TA, CD, MNND, and SDNND, from 16 different reef sites. Data file also contains three reef variables: LFTADen (Lionfish density /100m2), HASAve (averaged score of habitat complexity), and lionfish removal treatment (binary, yes/no).

``` r
d <- read.csv("data/FullCommNoLF.csv")
```

Quick model
===========

``` r
d_logged <- group_by(d, Site) %>%
  mutate_each(funs(log), dNr:SDNND)

d_means <- d_logged %>%
  summarise_each(funs(mean))

d_sds <- d_logged %>%
  summarise_each(funs(sd), dNr:SDNND) %>%
  dplyr::select(-Site)
names(d_sds) <- paste0(names(d_sds), "_sd")

d_sum <- data.frame(d_means, d_sds)

d_sum <- d_sum %>% mutate(HASAve = arm::rescale(HASAve), 
  LFTADen = arm::rescale(LFTADen),
  RemovTreat = arm::rescale(RemovTreat))

d_sum <- d_sum %>% left_join(readr::read_csv("data/SiteCoords.csv"))

# reshape2::melt(d, )

m1 <- lm(dNr ~ HASAve + LFTADen * RemovTreat, data = d_sum)
arm::display(m1)
#> lm(formula = dNr ~ HASAve + LFTADen * RemovTreat, data = d_sum)
#>                    coef.est coef.se
#> (Intercept)         1.27     0.02  
#> HASAve              0.14     0.05  
#> LFTADen            -0.03     0.05  
#> RemovTreat         -0.04     0.04  
#> LFTADen:RemovTreat  0.22     0.10  
#> ---
#> n = 16, k = 5
#> residual sd = 0.09, R-Squared = 0.60
```

Stan
====

Same model to check:

``` r
X <- model.matrix(dNr ~ 0 + HASAve + LFTADen * RemovTreat, data = d_sum)
stan_dat <- list(y_meas = d_sum$dNr, tau = d_sum$dNr_sd, N = nrow(d_sum), K = ncol(X), X = X)
```

``` r
writeLines(readLines("tnm.stan"))
#> data {
#>   int N;               // number of observations
#>   int K;               // number of predictors
#>   real y_meas[N];      // measurement of y
#>   matrix[N, K] X;      // model predictor matrix
#> }
#> parameters {
#>   vector[K] beta;      // vector of predictors
#>   real alpha;          // intercept
#>   real<lower=0> sigma; // residual sd
#> } 
#> model { 
#>   sigma ~ student_t(3, 0, 1);  // prior
#>   alpha ~ student_t(3, 0, 10);  // prior
#>   beta ~ student_t(3, 0, 1);   // prior
#>   y_meas ~ normal(alpha + X * beta, sigma); // likelihood
#> }
```

``` r
m_basic <- stan("tnm.stan", data = stan_dat, control = ctl)
```

``` r
m_basic
#> Inference for Stan model: tnm.
#> 4 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=4000.
#> 
#>          mean se_mean   sd  2.5%   25%   50%   75% 97.5% n_eff Rhat
#> beta[1]  0.14    0.00 0.06  0.04  0.11  0.14  0.18  0.25  2978    1
#> beta[2] -0.03    0.00 0.06 -0.15 -0.07 -0.03  0.01  0.09  3022    1
#> beta[3] -0.04    0.00 0.05 -0.15 -0.08 -0.04 -0.01  0.05  2916    1
#> beta[4]  0.22    0.00 0.11  0.00  0.15  0.22  0.29  0.44  2889    1
#> alpha    1.27    0.00 0.03  1.22  1.26  1.27  1.29  1.32  3000    1
#> sigma    0.10    0.00 0.02  0.07  0.08  0.10  0.11  0.16  1766    1
#> lp__    27.56    0.07 2.18 22.27 26.41 28.00 29.17 30.51  1109    1
#> 
#> Samples were drawn using NUTS(diag_e) at Sat Dec 10 18:16:03 2016.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

That looks the same to me.

Measurement error model
=======================

Now let's use the following model that allows for measurement error on y:

``` r
writeLines(readLines("tnm-meas.stan"))
#> data {
#>   int N;               // number of observations
#>   int K;               // number of predictors
#>   real y_meas[N];      // measurement of y
#>   real<lower=0> tau[N];   // measurement sd on y
#>   matrix[N, K] X;      // model predictor matrix
#> }
#> parameters {
#>   vector[K] beta;      // vector of predictors
#>   real alpha;          // intercept
#>   real<lower=0> sigma; // residual sd
#>   real y_raw[N];
#> }
#> transformed parameters {
#>   real y[N];           // unknown true y value
#>   for (i in 1:N) {
#>     y[i] = alpha + X[i, ] * beta + sigma * y_raw[i]; // non-centered parameterization 'trick'
#>   }
#> } 
#> model { 
#>   sigma ~ student_t(3, 0, 1);  // prior
#>   alpha ~ student_t(3, 0, 10);  // prior
#>   beta ~ student_t(3, 0, 1);   // prior
#>   y_meas ~ normal(y, tau);     // measurement model
#>   y_raw ~ normal(0, 1);        // non-centered parameterization 'trick'
#> }
#> generated quantities{
#>   real y_pred[N];
#>   for (i in 1:N) {
#>     y_pred[i] = alpha + X[i, ] * beta;
#>   }
#> }
```

``` r
m_meas <- stan("tnm-meas.stan", data = stan_dat, 
  pars = c("y", "y_raw"), include = FALSE, control = ctl)
```

``` r
m_meas
#> Inference for Stan model: tnm-meas.
#> 4 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=4000.
#> 
#>              mean se_mean   sd   2.5%    25%    50%    75% 97.5% n_eff
#> beta[1]      0.15    0.00 0.09  -0.03   0.09   0.15   0.21  0.33  4000
#> beta[2]     -0.06    0.00 0.13  -0.30  -0.15  -0.06   0.03  0.20  4000
#> beta[3]     -0.05    0.00 0.09  -0.22  -0.11  -0.05   0.00  0.12  4000
#> beta[4]      0.21    0.00 0.23  -0.25   0.05   0.21   0.36  0.65  4000
#> alpha        1.29    0.00 0.05   1.20   1.26   1.29   1.32  1.38  4000
#> sigma        0.05    0.00 0.04   0.00   0.02   0.04   0.07  0.16  2281
#> y_pred[1]    1.38    0.00 0.08   1.23   1.33   1.38   1.43  1.53  4000
#> y_pred[2]    1.28    0.00 0.13   1.04   1.20   1.28   1.36  1.53  4000
#> y_pred[3]    1.12    0.00 0.19   0.76   1.00   1.12   1.25  1.50  4000
#> y_pred[4]    1.36    0.00 0.09   1.18   1.30   1.36   1.42  1.55  4000
#> y_pred[5]    1.40    0.00 0.08   1.24   1.35   1.40   1.45  1.56  4000
#> y_pred[6]    1.36    0.00 0.09   1.19   1.30   1.36   1.42  1.53  4000
#> y_pred[7]    1.30    0.00 0.07   1.17   1.25   1.30   1.35  1.43  4000
#> y_pred[8]    1.33    0.00 0.06   1.22   1.29   1.33   1.37  1.44  4000
#> y_pred[9]    1.16    0.00 0.09   0.97   1.10   1.16   1.23  1.35  4000
#> y_pred[10]   1.35    0.00 0.06   1.22   1.31   1.35   1.39  1.47  4000
#> y_pred[11]   1.23    0.00 0.10   1.04   1.16   1.23   1.29  1.42  4000
#> y_pred[12]   1.12    0.00 0.13   0.87   1.04   1.12   1.21  1.37  4000
#> y_pred[13]   1.36    0.00 0.08   1.20   1.31   1.36   1.41  1.52  4000
#> y_pred[14]   1.23    0.00 0.07   1.09   1.18   1.23   1.28  1.38  4000
#> y_pred[15]   1.19    0.00 0.10   1.00   1.13   1.19   1.26  1.38  4000
#> y_pred[16]   1.37    0.00 0.17   1.03   1.25   1.37   1.48  1.70  4000
#> lp__       -15.33    0.09 3.38 -22.71 -17.49 -15.09 -12.95 -9.50  1326
#>            Rhat
#> beta[1]       1
#> beta[2]       1
#> beta[3]       1
#> beta[4]       1
#> alpha         1
#> sigma         1
#> y_pred[1]     1
#> y_pred[2]     1
#> y_pred[3]     1
#> y_pred[4]     1
#> y_pred[5]     1
#> y_pred[6]     1
#> y_pred[7]     1
#> y_pred[8]     1
#> y_pred[9]     1
#> y_pred[10]    1
#> y_pred[11]    1
#> y_pred[12]    1
#> y_pred[13]    1
#> y_pred[14]    1
#> y_pred[15]    1
#> y_pred[16]    1
#> lp__          1
#> 
#> Samples were drawn using NUTS(diag_e) at Sat Dec 10 18:16:08 2016.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

We can inspect and plot the output:

``` r
posterior <- extract(m_meas, inc_warmup = FALSE, permuted = FALSE, pars = "y_pred", 
  include = FALSE)
mcmc_trace(posterior)
```

![](tnm_files/figure-markdown_github/stan-check-plots-1.png)

``` r
names(as.data.frame(X))
#> [1] "HASAve"             "LFTADen"            "RemovTreat"        
#> [4] "LFTADen:RemovTreat"

mcmc_areas(as.matrix(m_meas), regex_pars = "beta")
```

![](tnm_files/figure-markdown_github/stan-check-plots-2.png)

``` r
mcmc_areas(as.matrix(m_basic), regex_pars = "beta")
```

![](tnm_files/figure-markdown_github/stan-check-plots-3.png)

So the posteriors are definitely wider in the case when we allow for the measurement error.

Other responses
===============

Now that the above models are working, let's apply them to the various responses.

The following function will format the data and fit the model for a given response. It also has the option of centering or not centering the removal variable.

``` r
fit_tnm <- function(response) {
  f <- paste(response, "~ 0 + HASAve + LFTADen * RemovTreat")
  X <- model.matrix(as.formula(f), data = d_sum)
  stan_dat <- list(y_meas = d_sum[, response], tau = d_sum$dNr_sd, N = nrow(d_sum), 
    K = ncol(X), X = X)
  m <- stan("tnm-meas.stan", data = stan_dat,
    pars = c("y", "y_raw"), include = FALSE, control = list(adapt_delta = 0.99), 
    iter = 3000, chains = 4)
  m
}
```

Now let's apply the function to each of the responses.

``` r
responses <- names(d)[2:7]
out_cent <- lapply(responses, fit_tnm)
names(out_cent) <- responses
```

Here we will extract the posteriors from the models and reformat the output for plotting. Also, we will calculate the effect of lionfish density for the case without removals (-0.5) and for the case with removals (0.5).

``` r
stopifnot(identical(mean(d$RemovTreat), 0.5)) # Just in case! 
p <- plyr::ldply(out_cent, stanhelpers::extract_df, output = "wide_df")
pred <- select(p, .id, starts_with("y_pred")) %>% rename(response = `.id`)
est <- select(p, .id, starts_with("beta")) %>% rename(response = `.id`)
names(est)[2:5] <- names(as.data.frame(X))
est <- mutate(est, LFTADen_w_removal = LFTADen + 0.5 * `LFTADen:RemovTreat`,
  LFTADen_no_removal = LFTADen - 0.5 * `LFTADen:RemovTreat`)
est <- reshape2::melt(est) # make a long format for ggplot
```

``` r
labs <- c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5)
ggplot(est, aes(variable, value, fill = response, colour = response)) + 
  geom_hline(yintercept = 0, lty = 2) + xlab("") +
  scale_y_continuous(breaks = log(labs), labels = labs, limits = range(log(labs))) +
  geom_violin(position = position_dodge(width = 0.8), alpha = 0.5,
    draw_quantiles = 0.5) +
  coord_flip() +
  viridis::scale_fill_viridis(discrete = TRUE) +
  viridis::scale_color_viridis(discrete = TRUE) +
  theme_light() +
  ylab("Coefficient estimate")
```

![](tnm_files/figure-markdown_github/coefficient-plot-1.png)

``` r
ggsave("figs/tnm-estimates.pdf", width = 6, height = 6.5)
```

In the above plot, I labeled the x axis with the exponentiated versions of the coefficients. These are the multiplicative effects. So, for example, a value of 0.6 means that the response will be 60% of what it was if the predictor increases by 2 standard deviations (or in the case of the treatment, the lionfish are removed).

We can calculate the probability a given coefficient is above or less than 0 (i.e. the multiplicative effect is above or below 1) (because this is a Bayesian model). You can use these values when you report results. We can also calculate lots of other things depending on what would be meaningful (e.g. credible intervals on the effects).

``` r
sum_table <- est %>% group_by(variable, response) %>% 
  summarize(
    prob_less_0 = round(sum(value < 0)/n(), 2),
    prob_above_0 = round(sum(value > 0)/n(), 2))
knitr::kable(sum_table)
```

| variable             | response |  prob\_less\_0|  prob\_above\_0|
|:---------------------|:---------|--------------:|---------------:|
| HASAve               | CD       |           0.25|            0.75|
| HASAve               | dCr      |           0.62|            0.38|
| HASAve               | dNr      |           0.06|            0.94|
| HASAve               | MNND     |           0.01|            0.99|
| HASAve               | SDNND    |           0.54|            0.46|
| HASAve               | TA       |           0.30|            0.70|
| LFTADen              | CD       |           0.73|            0.27|
| LFTADen              | dCr      |           0.50|            0.50|
| LFTADen              | dNr      |           0.67|            0.33|
| LFTADen              | MNND     |           0.87|            0.13|
| LFTADen              | SDNND    |           0.39|            0.61|
| LFTADen              | TA       |           0.66|            0.34|
| RemovTreat           | CD       |           0.87|            0.13|
| RemovTreat           | dCr      |           0.88|            0.12|
| RemovTreat           | dNr      |           0.75|            0.25|
| RemovTreat           | MNND     |           0.66|            0.34|
| RemovTreat           | SDNND    |           0.97|            0.03|
| RemovTreat           | TA       |           0.95|            0.05|
| LFTADen:RemovTreat   | CD       |           0.47|            0.53|
| LFTADen:RemovTreat   | dCr      |           0.40|            0.60|
| LFTADen:RemovTreat   | dNr      |           0.19|            0.81|
| LFTADen:RemovTreat   | MNND     |           0.47|            0.53|
| LFTADen:RemovTreat   | SDNND    |           0.07|            0.93|
| LFTADen:RemovTreat   | TA       |           0.09|            0.91|
| LFTADen\_w\_removal  | CD       |           0.64|            0.36|
| LFTADen\_w\_removal  | dCr      |           0.45|            0.55|
| LFTADen\_w\_removal  | dNr      |           0.41|            0.59|
| LFTADen\_w\_removal  | MNND     |           0.75|            0.25|
| LFTADen\_w\_removal  | SDNND    |           0.16|            0.84|
| LFTADen\_w\_removal  | TA       |           0.32|            0.68|
| LFTADen\_no\_removal | CD       |           0.73|            0.27|
| LFTADen\_no\_removal | dCr      |           0.58|            0.42|
| LFTADen\_no\_removal | dNr      |           0.87|            0.13|
| LFTADen\_no\_removal | MNND     |           0.86|            0.14|
| LFTADen\_no\_removal | SDNND    |           0.84|            0.16|
| LFTADen\_no\_removal | TA       |           0.93|            0.07|

Here's what I see: *(numbers might have changed slightly)*

-   None of these effects are overly strong
-   There's a reasonably high probability that most of the responses are lower in the case of lionfish removals (for a site with average lionfish density) (see the `RemovTreat` effect). For example, there is about a 95% probability this is true for `TA` and about a 88% probability this is true for `CD` and `dCr`.
-   There is a fairly high probability (ranging from 0.6 to 0.93) that some of these responses (TA, dCr, dNr) are negatively related with lionfish density (`LFTADen`). This is for an average site, or in other words across all sites including those with and without removals.
-   The effect of lionfish density looks a bit stronger (negative) in the case of no removals `LFTADen_no_removal`, BUT the interaction between lionfish density and treatment is very weak (or at least very uncertain) (see `LFTADen:RemovTreat`).
-   There is weak evidence for an effect of the habitat variable on the responses with the exception of positive relationship between the habitat variables and `dNr` and `MMND`, with ~0.98 or 0.99 probability for the latter. (`HASAve`)

Residuals
=========

``` r
p <- reshape2::melt(pred, variable.name = "site_i", value.name = "predicted")
p <- p %>% group_by(response, site_i) %>%
  summarise(predicted = mean(predicted))

obs <- select(d_sum, Site:SDNND, lat, lon) %>%
  reshape2::melt(id.vars = c("Site", "lat", "lon"), variable.name = "response", 
    value.name = "observed")

stopifnot(identical(nrow(d_sum), length(unique(p$site_i)))) # just to check
lookup <- data_frame(Site = d_sum$Site, site_i = paste0("y_pred_", seq_len(nrow(d_sum))))
res <- inner_join(p, lookup) %>% inner_join(obs) %>%
  mutate(residual = observed - predicted)

ggplot(res, aes(predicted, observed)) + geom_point() +
  facet_wrap(~response, scales = "free") +
  geom_abline(intercept = 0, slope = 1, lty = 2)
```

![](tnm_files/figure-markdown_github/residuals-1.png)

``` r

ggplot(res, aes(predicted, residual)) + geom_point() +
  facet_wrap(~response, scales = "free_x") +
  geom_abline(intercept = 0, slope = 0, lty = 2)
```

![](tnm_files/figure-markdown_github/residuals-2.png)

``` r
ggsave("figs/fitted-vs-obs-tnm-residuals.pdf", width = 9, height = 6.5)

ggplot(res, aes(lon, lat, colour = residual)) + geom_point(size = 4) +
  facet_wrap(~response) +
  scale_color_gradient2()
```

![](tnm_files/figure-markdown_github/residuals-3.png)

``` r
ggsave("figs/spatial-tnm-residuals.pdf", width = 10, height = 6.5)
```

Looks pretty good.

Check mixing
============

``` r
color_scheme_set("mix-blue-red")
lapply(out_cent, function(x) mcmc_trace(extract(x, inc_warmup = FALSE, permuted = FALSE, 
  pars = c("y_pred"), include = FALSE)))
#> $dNr
```

![](tnm_files/figure-markdown_github/unnamed-chunk-16-1.png)

    #> 
    #> $dCr

![](tnm_files/figure-markdown_github/unnamed-chunk-16-2.png)

    #> 
    #> $TA

![](tnm_files/figure-markdown_github/unnamed-chunk-16-3.png)

    #> 
    #> $CD

![](tnm_files/figure-markdown_github/unnamed-chunk-16-4.png)

    #> 
    #> $MNND

![](tnm_files/figure-markdown_github/unnamed-chunk-16-5.png)

    #> 
    #> $SDNND

![](tnm_files/figure-markdown_github/unnamed-chunk-16-6.png)

``` r

lapply(out_cent, function(x) {
  broom::tidyMCMC(x, rhat = TRUE, ess = TRUE) %>% 
    filter(!grepl("y", term)) %>%
    select(term, rhat, ess) %>%
    mutate(rhat = round(rhat, 2), ess = round(ess, 2))
})
#> $dNr
#>      term rhat  ess
#> 1 beta[1]    1 6000
#> 2 beta[2]    1 6000
#> 3 beta[3]    1 6000
#> 4 beta[4]    1 6000
#> 5   alpha    1 6000
#> 6   sigma    1 3128
#> 
#> $dCr
#>      term rhat  ess
#> 1 beta[1]    1 6000
#> 2 beta[2]    1 6000
#> 3 beta[3]    1 6000
#> 4 beta[4]    1 6000
#> 5   alpha    1 6000
#> 6   sigma    1 2588
#> 
#> $TA
#>      term rhat  ess
#> 1 beta[1]    1 3902
#> 2 beta[2]    1 4040
#> 3 beta[3]    1 4229
#> 4 beta[4]    1 4651
#> 5   alpha    1 4410
#> 6   sigma    1 1900
#> 
#> $CD
#>      term rhat  ess
#> 1 beta[1]    1 5329
#> 2 beta[2]    1 6000
#> 3 beta[3]    1 6000
#> 4 beta[4]    1 6000
#> 5   alpha    1 6000
#> 6   sigma    1 2516
#> 
#> $MNND
#>      term rhat  ess
#> 1 beta[1]    1 2863
#> 2 beta[2]    1 3600
#> 3 beta[3]    1 3171
#> 4 beta[4]    1 4383
#> 5   alpha    1 3244
#> 6   sigma    1 1864
#> 
#> $SDNND
#>      term rhat  ess
#> 1 beta[1]    1 3008
#> 2 beta[2]    1 3578
#> 3 beta[3]    1 3099
#> 4 beta[4]    1 3944
#> 5   alpha    1 3212
#> 6   sigma    1 1656
```

Looks great.

brms package checks
===================

Similar, but subtly different models as checks (ignore).

``` r
# library(brms)
# m2 <- brm(dNr ~ 1 + HASAve + LFTADen * RemovTreat, data = d_sum)

# m3 <- brm(dNr | se(dNr_sd) ~ 1 + HASAve + LFTADen * RemovTreat, data = d_sum,
#   prior = c(prior(student_t(5, 0, 2), class = "b"), prior(student_t(5, 0, 5), class = "Intercept")))

# m4 <- brm(dNr | se(dNr_sd) ~ 1 + HASAve + LFTADen * RemovTreat + (1|Site), data = d_sum,
#   control = list(adapt_delta = 0.99),
#   prior = c(prior(student_t(5, 0, 2), class = "b"),
#     prior(student_t(5, 0, 5), class = "Intercept"),
#     prior(student_t(5, 0, 2), class = "sd")))
```
