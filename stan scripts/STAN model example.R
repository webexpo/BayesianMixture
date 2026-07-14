############################################################################
#
#  Example to run the STAN models for the bernouilli lognormal  and censored models
#  
#
#
##############################################################################

###############


library(LaplacesDemon)
library(rstan)
library(bayesplot)
library(ggplot2)
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
stanmodel.trunc


###############  DATA GENERATION 

    data0 <- fun.create.data.singleLOD(   prop.true.zero=0.3 ,   # proportion of true zeroes (omega)
                                          n=100,                  # total sample size
                                          gm=100,                 # geometric mean of the lognormal distribution part of mixture
                                          gsd=2.5,                # geometric standard deviation of the lognormal distribution part of mixture
                                          prop.nonzero.censored=0.01  # proportion of the lognormal distribution supposed to be censored
    )


###############################  generating posteriors
    
    
    post.trunc <- fun.stan.trunc( data_sample = data0 , n.iter = 5000 , n.warmup = 2500 , n.chain = 4 ) 
      
    post.mixt <- fun.stan.mix( data_sample = data0 , n.iter = 5000 , n.warmup = 2500 , n.chain = 4 )  

    post.trunc.2 <- fun.stan.cens( data_sample = data0 , n.iter = 5000 , n.warmup = 2500 , n.chain = 4 )  
    
  
#### bayesian analysis


mixture.analysis <- stan_bayes.mixt( data.sample = data0 , n.iter = 4000 , n.warmup = 2000 , mymodel=stanmodel.mixt)
  
object.size(mixture.analysis)

censored.analysis <- stan_bayes.cens( data.sample = data0 , n.iter = 4000 , n.warmup = 2000  , mymodel=stanmodel.cens)

object.size(censored.analysis)


####  point estimates and confidence intervals in a table


tab.res <- estimates( mixture.analysis = mixture.analysis, censored.analysis =  censored.analysis , conf = 89)

tab.res <- othermetrics( mixture.analysis = mixture.analysis, censored.analysis =  censored.analysis , conf = 89 , oel=100 , target_perc = 95/100)





####  Odds of the mixture model being better than the censored model

GOF <- goodness.of.fit( mixture.analysis = mixture.analysis, censored.analysis =  censored.analysis )


###### the MCMC chains

chains <- posterior.chains( mixture.analysis = mixture.analysis, censored.analysis =  censored.analysis  ) 
  


#### posterior predictive comparison of the censored and mixture models


pp.graph <- posterior.pred.graphs( post.results = chains , nrep = 5 , data0 = data0 ) 

pp.graph$p.mix
pp.graph$p.cens


#################################################### posterion predictions - quantties

pp.quant <- posterior.pred.quant( post.results  = chains, 
                                  nrep = 500 , data0 = data0 , 
                                  type = "propcens" )
pp.quant$p.mix  
pp.quant$p.cens  


pp.quant <- posterior.pred.quant( post.results  = chains, 
                                  nrep = 500 , data0 = data0 , 
                                  type = "p95" )
pp.quant$p.mix  
pp.quant$p.cens  
