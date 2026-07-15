# BayesianMixture

R scripts accompanying the paper:

> Lavoue J, Burstyn I. (2021) Evidence of Absence: Bayesian Way to Reveal True Zeros Among Occupational Exposures. *Annals of Work Exposures and Health*, 65(1), 84–95. https://doi.org/10.1093/annweh/wxaa086

## Contents

This repository contains the code used to estimate the zero-inflated lognormal mixture model (proportion of true zeros, geometric mean, and geometric standard deviation) presented in the paper, including:
- the truncated model
- the censored Bernoulli-lognormal mixture model
- the standard censored lognormal model (baseline)

## Implementation

Bayesian estimation (MCMC) implemented in both **JAGS** and **STAN**.

An online application allowing reproduction of the analyses is available here:
https://lavoue.shinyapps.io/Mixture_app_STAN/

Note : as of july 2026, not all calculus implemented in the Web app are available here

## Citation

If you use this code, please cite the paper above.
