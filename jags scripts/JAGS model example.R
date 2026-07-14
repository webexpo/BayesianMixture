############################################################################
#
#  Example to run the JAGS models for the bernouilli lognormal  and truncated mixture
#  
#
#
##############################################################################

###############

library(rjags)
library(bayesplot)
library(coda)
library(rethinking)
library(LaplacesDemon)


############# file sourcing
setwd("F:/Dropbox/bureau/RStudio/BayesianMixtureModels/Mixture2019")

setwd("D:/Dropbox/bureau/RStudio/BayesianMixtureModels/Mixture2019")


setwd("C:/jerome/Dropbox/bureau/RStudio/BayesianMixtureModels/Mixture2019")

source("Ancillary functions/data generation script.R")

source("JAGS related functions/JAGS final trunc and mixture model statements.R")

source("JAGS related functions/JAGS model functions.R")

###############  DATA GENERATION 

    data0 <- fun.create.data.singleLOD(   prop.true.zero=0.25 ,   # proportion of true zeroes (omega)
                                          n=30,                  # total sample size
                                          gm=300,                 # geometric mean of the lognormal distribution part of mixture
                                          gsd=2.5,                # geometric standard deviation of the lognormal distribution part of mixture
                                          prop.nonzero.censored=0.25  # proportion of the lognormal distribution supposed to be censored
    )


###############################  generating posteriors
    
    
    post.trunc <- fun.jags( data.sample = data0 , n.iter = 50000 , n.burnin = 5000 , type = "truncated" ) 
      
    post.mixt <- fun.jags( data.sample = data0 , n.iter = 50000 , n.burnin = 5000 , type = "mixture" ) 


###############################  comparing the poinst estimates and CIs
    
    
    fun.comp.trunc.mixt( mixt.results = post.mixt ,
                         trunc.results = post.trunc ,
                         data0 = data0)
    
    
###############################  comparing the densities   
    
       fun.compare.dens(chain.1 = post.mixt$gm ,
                        chain.2 = post.trunc$gm ,
                        true.value = data0$trueGM ,
                        parameter = "GM")
      
       fun.compare.dens(chain.1 = post.mixt$gsd ,
                        chain.2 = post.trunc$gsd ,
                        true.value = data0$trueGSD ,
                        parameter = "GSD")
    
       fun.compare.dens(chain.1 = post.mixt$omega ,
                        chain.2 = post.trunc$omega ,
                        true.value = data0$omega ,
                        parameter = "omega")
    
    
######## diagnostics with JAGS MIXTURE
       
       ####### PREPARING DATA FOR JAGS
       
       # total sample size
       
       n <- data0$n
       
       # number of unobserved values
       
       n0 <- data0$n0
       
       #constant for the ones tricks procedure in JAGS, which allows writing the log-likelihood in full
       
       C <-100000
       
       zeros <- rep( 0 , n )
       
       # Pi conrtant
       
       pi.constant <-3.14159265359
       
       #initial values
       
       inits <- list( mu = log(0.3) , log.sigma = log(2.5) , omega = 0 ) # arbitrary but do not matter in the calculation
       # alternatively one can use random numbers from the priors
       #list containing the data for JAGS
       
       jags.data <- list(x = data0$x,
                         is.observed = data0$x.is.observed,
                         n = n,
                         n0 = n0,
                         C = C,
                         zeros = zeros,
                         pi.constant = pi.constant )
       
       ###### GENERATING POSTERIOR SAMPLES
       
       # model initialization
       
       j.mod <- jags.model( file = textConnection( model.jags.mixture ),
                            data = jags.data,
                            inits = inits,
                            n.chains = 2,
                            n.adapt = 100)
       
       #burnin iterations
       
       update(j.mod,n.iter = 5000)
       
       ###### generating the MCMC chains
       
       c.out <-coda.samples( model = j.mod,
                             variable.names = c('mu','sigma','omega',"ll"),   #,"ll"),
                             n.iter = 10000,
                             thin=1)
       j.out <-jags.samples( model = j.mod,
                             variable.names = c('mu','sigma','omega'),
                             n.iter = 50000,
                             thin=1)
       
       
       
    
        ######### extracting the posterior samples for compatibility with bayesplot
        
        j.array.mixt <- array( dim = c( 50000 , 2 , 3) )
        
        j.array.mixt[,,1] <-j.out$mu[1,,]
        j.array.mixt[,,2] <-j.out$sigma[1,,]
        j.array.mixt[,,3] <-j.out$omega[1,,]
        
        dimnames(j.array.mixt)[[3]] <- c("mu", "sigma", "omega")
        
        mcmc_combo(j.array.mixt, combo = c("dens", "trace"))
        
        mcmc_acf_bar(j.array.mixt)
        
        gelman.diag(c.out[ , c("mu","sigma","omega") ] )
        raftery.diag(c.out[ , c("mu","sigma","omega") ])
        gelman.plot(c.out)
        autocorr(c.out)
        crosscorr.plot(c.out[ , c("mu","sigma","omega") ])
        cumuplot(c.out)
        effectiveSize(c.out[ , c("mu","sigma","omega") ])
        geweke.diag(c.out[ , c("mu","sigma","omega") ])
        heidel.diag(c.out[ , c("mu","sigma","omega") ])
        
        
        
        loglik <- matrix( nrow = 50000  , ncol = data0$n )
        
        for (i in 1:data0$n ) loglik[,i] <-c.out[[2]][,paste("ll[",i,"]",sep="")]
        
        WAIC( x = loglik )
        
        



