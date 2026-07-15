######################################
#
#  generating data from the bernouilli lognormal mixture
#  
#  V0.6 July / 2026
#
######################################


########data creation function

fun.create.data.singleLOD <- function(prop.true.zero=0.15 ,   # proportion of true zeroes (omega)
                                      n=200,                  # total sample size
                                      gm=100,                 # geometric mean of the lognormal distribution part of mixture
                                      gsd=2.5,                # geometric standard deviation of the lognormal distribution part of mixture
                                      prop.nonzero.censored=0.15  # proportion of the lognormal distribution supposed to be censored
)  {
  
  
  #Limit of quantification as a function of gm, gsd and prop.nonzero.censored 
  
  LOQ <-exp( qnorm( prop.nonzero.censored  , log(gm) , log(gsd) ) )
  
  ## non zero data generation (generation of lognormal data) for the non zero proportion of the sample
  
  omega.sample <-  rbeta(1 , trunc(n*prop.true.zero) , trunc(n*(1-prop.true.zero)) )
  
  n.non.zero <- round( n * (1 - omega.sample))
  
  data0 <- rlnorm( n = n.non.zero , mean = log( gm ), sd = log( gsd ) )
  
  data0.observed.status <- rep( TRUE , length(data0) )
  
  ## sample statistics
  
  true_sampleGM <- exp(mean(log(data0)))
  
  true_sampleGSD <- exp(sd(log(data0)))
  
  true_propcens <- length(data0[data0<LOQ])/length(data0)
  
  #censoring the non zero observations
  
  data0.observed.status[data0<LOQ] <- FALSE
  
  #adding the true zeroes
  
  data0<-c( data0 , rep( 0 , n - length( data0 ) ) )
  
  data0.observed.status <- c( data0.observed.status , rep( FALSE , n - length( data0.observed.status ) ) )
  
  # n1 in Taylor et al.
  
  n1 <- length( data0.observed.status[ data0.observed.status ] )
  
  # n0 in Taylor et al.
  
  n0 <- n - n1
  
  # replacing unobserved values with LOQ
  
  data0[!data0.observed.status] <- LOQ
  
  # shuffling the data
  
  index <-sample(1:n)
  
  data0 <-data0[index]
  
  data0.observed.status <-  data0.observed.status[index]
  
  ## in the expostat format
  
  expostats.x <- as.character(signif(data0,3))
  expostats.x[ !data0.observed.status ] <- paste( "<" , expostats.x[!data0.observed.status] , sep = "")
  
  
  ### final data
  
  return(list( x=data0,                                    # final vector of values, fixed at censoring point when unobserved
               x.is.observed=data0.observed.status,        # logical vector, TRUE for observed records
               true_sampleGM=true_sampleGM,                # GM of the original lognormal sample
               true_sampleGSD=true_sampleGSD,              # GSD of the original lognormal sample
               trueGM=gm,                                  # true GM
               trueGSD=gsd,                                 # true GSD               
               true_propcens=true_propcens,                # censorship proportion in the original lognormal sample
               omega = prop.true.zero,                     # selected proportion of true zeroes
               omega.sample = omega.sample,                 # # proportion of true zeroes in the sample
               n = n,                                      # total sample size
               n1 = n1,                                    # number of observed values
               n0 = n0,                                    # number of unobserved values
               expostats.x = expostats.x,                   #expostats format
               LOQ = LOQ))                                  #Censoring point 
  
}
