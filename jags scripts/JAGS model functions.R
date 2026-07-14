#######################################################################################
#
#
#  Function that creates posteriors for the TRUNCATED and MIXTURE models
#
#
#######################################################################################

###### inputs

  # data.sample : data formatted in the same form as output by the fun.create.data.singleLOD function

  # n.iter number of iterations

  # n.burnin burnin iterations

  # type : "truncated" or "mixture"

###### outputs

  # chain for gm

  # chain for gsd

  # chain for omega


###### required scripts :

  #  JAGS final trunc and mixture model statements.R


fun.jags <- function( data.sample , n.iter = 25000 , n.burnin = 5000 , type = "mixture" ) 
  
        {
  
  ####### fitting the mixture model
  if ( type == "mixture" )
          
                  {
                   
          ####### PREPARING DATA FOR JAGS
          
          # total sample size
          
          n <- data.sample$n
          
          # number of unobserved values
          
          n0 <- data.sample$n0
          
          #constant for the ones tricks procedure in JAGS, which allows writing the log-likelihood in full
          
          C <-100000
          
          zeros <- rep( 0 , n )
          
          # Pi conrtant
          
          pi.constant <-3.14159265359
          
          #initial values
          
          inits <- list( mu = log(0.3) , log.sigma = log(2.5) , omega = 0 ) # arbitrary but do not matter in the calculation
          # alternatively one can use random numbers from the priors
          #list containing the data for JAGS
          
          jags.data <- list(x = data.sample$x,
                            is.observed = data.sample$x.is.observed,
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
                               n.chains = 1,
                               n.adapt = 100)
          
          #burnin iterations
          
          update(j.mod,n.iter = n.burnin)
          
          ###### generating the MCMC chains
          
          j.out <-coda.samples( model = j.mod,
                                variable.names = c('mu','sigma','omega'),
                                n.iter = n.iter,
                                thin=1)
          
          ######### extracting the posterior samples
          
          gm.chain <- as.numeric( exp( j.out[[1]][ , "mu" ] ) )
          gsd.chain <- as.numeric( exp( j.out[[1]][ , "sigma" ]) )
          omega.chain <- as.numeric( j.out[[1]][ , "omega" ] )
           
                          
  }
  
  ####### fitting the truncated model
  if ( type == "truncated" )
    
  {
  
          ####### PREPARING DATA FOR JAGS      
          
          
          # total sample size
          
          n <- length(data.sample$x)
          
          # number of unobserved values
          
          n0 <- data.sample$n0
          
          # truncation point : Xd
          
          Xd <- data.sample$LOQ
          
          #list containing the data for JAGS
          
          jags.data = list( X1 = data.sample$x[ data.sample$x.is.observed ], # observed values above LOD
                            n0 = n0,                           # number of unobserved values
                            n = n,                             # total sample size
                            Xd = Xd )                          # truncation point
          
          #initial values
          
          inits <- list( mu = log(0.3) , log.sigma = log(2.5) ) # arbitrary but do not matter in the calculation
          # alternatively one can use random numbers from the priors
          
          
          ###### GENERATING POSTERIOR SAMPLES
          
          # model initialization
          
          j.mod <- jags.model(file = textConnection(model.jags.truncated),
                              data =  jags.data,
                              inits = inits,
                              n.chains = 1,
                              n.adapt = 100)
          
          #burnin iterations
          
          update(j.mod,n.iter = n.burnin)
          
          
          ###### generating the MCMC chains
          
          j.out <-coda.samples( model = j.mod,
                                variable.names = c('mu','sigma','omega'),
                                n.iter = n.iter,
                                thin=1)
          
          ######### extracting the posterior samples
          
          #unconstrained
          gm.chain.unc <- as.numeric( exp( j.out[[1]][ , "mu" ] ) )
          gsd.chain.unc <- as.numeric(exp( j.out[[1]][ , "sigma" ]) )
          omega.chain.unc <- as.numeric( j.out[[1]][ , "omega" ] )
          
          posterior <- data.frame( gm = gm.chain.unc , 
                                   gsd = gsd.chain.unc , 
                                   omega = omega.chain.unc )
          
          #restrict to omega>0 due to constraint on non-negative probability
          #restrict to omega =<n0/n due to constraint of observed data below Xd
          
          posterior <- subset( posterior , posterior$omega>=0 & posterior$omega<=n0/n )
          
          gm.chain <- posterior$gm
          gsd.chain <- posterior$gsd
          omega.chain <- posterior$omega 
          
    
  }
  
  result <- data.frame( gm = gm.chain ,
                        gsd = gsd.chain ,
                        omega = omega.chain)
  
  return( result )
          
        }
