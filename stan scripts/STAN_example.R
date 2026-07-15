############################################################################
#
#  Example to run the STAN models for the bernouilli lognormal  and censored models
#  
#  V0.6 July 2026
#
##############################################################################

# 1. Libraries and functions --------------------------------------------------------


library(LaplacesDemon)
library(rstan)
#library(bayesplot)
#library(ggplot2)
library(coda)
library(loo)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

############# file sourcing

source("stan scripts/data_generation.R")
source("stan scripts/data_preparation.R")
source("stan scripts/STAN functions.R")


############ the script assumes the following models have been compiled using the STAN_models.R script

stanmodel.cens 
stanmodel.mixt



# 2. DATA GENERATION -------------------------------------------------------- 

    data0 <- fun.create.data.singleLOD(   prop.true.zero=0.3 ,   # proportion of true zeroes (omega)
                                          n=100,                  # total sample size
                                          gm=100,                 # geometric mean of the lognormal distribution part of mixture
                                          gsd=2.5,                # geometric standard deviation of the lognormal distribution part of mixture
                                          prop.nonzero.censored=0.01  # proportion of the lognormal distribution supposed to be censored
    )


    
  
# 3. Bayesian analysis --------------------------------------------------------


mixture.analysis <- stan_bayes.mixt( data.sample = data0 , n.iter = 4000 , n.warmup = 2000 , mymodel=stanmodel.mixt)
  
object.size(mixture.analysis)

censored.analysis <- stan_bayes.cens( data.sample = data0 , n.iter = 4000 , n.warmup = 2000  , mymodel=stanmodel.cens)

object.size(censored.analysis)


# 4. Point estimates and confidence intervals in a table -----------------------------------------------


tab.res <- estimates( mixture.analysis = mixture.analysis, censored.analysis =  censored.analysis , conf = 90)

tab.res <- othermetrics( mixture.analysis = mixture.analysis, censored.analysis =  censored.analysis , conf = 90 , oel=100 , target_perc = 95/100)





# 5.  Odds of the mixture model being better than the censored model -----------------------------------------------

GOF <- goodness.of.fit( mixture.analysis = mixture.analysis, censored.analysis =  censored.analysis )


# 6. the MCMC chains ------------------------------------------------

chains <- posterior.chains( mixture.analysis = mixture.analysis, censored.analysis =  censored.analysis  ) 
  


