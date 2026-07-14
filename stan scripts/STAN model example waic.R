############################################################################
#
#  Example to compare the MIXT and CENSORED models
#  
#
##############################################################################


###############

library(rstan)
library(bayesplot)
library(coda)
library(rethinking)
library(loo)


############# file sourcing


setwd("C:/jerome/Dropbox/bureau/RStudio/BayesianMixtureModels/Mixture2019")

source("Ancillary functions/data generation script.R")

source("STAN related functions/STAN model functions.R")

source("STAN related functions/censoredVSmixture comparison functions.R")

stanmodel.cens <-readRDS("STAN related functions/stanmodel.cens.RDS")

stanmodel.mixt <-readRDS("STAN related functions/stanmodel.mixt.RDS")


###############  DATA GENERATION 

data0 <- fun.create.data.singleLOD(   prop.true.zero=0.1 ,   # proportion of true zeroes (omega)
                                      n=200,                  # total sample size
                                      gm=50,                 # geometric mean of the lognormal distribution part of mixture
                                      gsd=2,                # geometric standard deviation of the lognormal distribution part of mixture
                                      prop.nonzero.censored=0.3  # proportion of the lognormal distribution supposed to be censored
)


############## comparison


fun.compare.waic(data0)


