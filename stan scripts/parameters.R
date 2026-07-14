####################################
#
#
#  parameter estimates
#
#
########################


###### all models  : estimation of GM of the lognormal

gm <- function( mu.chain  , conf) {
  
  chain <- exp(mu.chain)
  
  chain <- exp(mu.chain)
  
  est <-median(chain)
  
  lcl <-as.numeric(HPDinterval( as.mcmc( as.numeric(chain )), prob = conf/100))[1]
  
  ucl <-as.numeric(HPDinterval( as.mcmc( as.numeric(chain ) ), prob = conf/100))[2]
  
  return(list(est=est,lcl=lcl,ucl=ucl))
  
  
}

###### all models  : estimation of GSD of the lognormal

gsd <- function( sigma.chain  , conf) {
  
  chain <- exp(sigma.chain)
  
  est <-median(chain)
  
  lcl <-as.numeric(HPDinterval( as.mcmc( as.numeric(chain )), prob = conf/100))[1]
  
  ucl <-as.numeric(HPDinterval( as.mcmc( as.numeric(chain ) ), prob = conf/100))[2]
  
  return(list(est=est,lcl=lcl,ucl=ucl))
  
  
}


###### mixture and trunc models  : estimation of omega 

omega.mixt <- function( omega.chain  , conf) {
  
  chain <- omega.chain
  
  est <-median(chain)
  
  lcl <-as.numeric(HPDinterval( as.mcmc( as.numeric(chain )), prob = conf/100))[1]
  
  ucl <-as.numeric(HPDinterval( as.mcmc( as.numeric(chain ) ), prob = conf/100))[2]
  
  return(list(est=est,lcl=lcl,ucl=ucl))
  
  
}


###### mixture and trunc  models  : estimation of the "target_perc" percentile 

perc.mixt <- function( mu.chain , sigma.chain , omega.chain , target_perc = 0.5 , conf) {
  

  quant <- ( target_perc - omega.chain) / ( 1 - omega.chain)
  
  chain <- numeric(length(mu.chain))
  
  criterion <- omega.chain >= target_perc
  
  chain[ criterion ] <- 0
  
  chain[ !criterion ] <- exp( qnorm( quant[ !criterion  ] , mu.chain[ !criterion  ] , sigma.chain[ !criterion  ] ) )
  
  est <-median(chain)
  
  lcl <-as.numeric(HPDinterval( as.mcmc( as.numeric(chain )), prob = conf/100))[1]
  
  ucl <-as.numeric(HPDinterval( as.mcmc( as.numeric(chain ) ), prob = conf/100))[2]
  
  return(list(est=est,lcl=lcl,ucl=ucl))
  
  
}


perc.mixt.true <- function( mu , sigma , omega , target_perc = 0.5 ) {
  
  
  if (omega >= target_perc) res <- 0
  
  if (omega < target_perc) {
    
    quant <- ( target_perc - omega) / ( 1 - omega)
    
    res <- exp( qnorm( quant , mu , sigma ) )
  }
  
  return(res)
  
  
}

###### mixture and trunc  models  : estimation of the exceeedance 

frac.mixt <- function( mu.chain , sigma.chain , omega.chain ,  conf , c.oel) {
  
  
  chain <-100*(1-pnorm((log((c.oel))-mu.chain)/sigma.chain))
  
  chain <- chain*(1-omega.chain)
  
  est <-median(chain)
  
  lcl <-as.numeric(HPDinterval( as.mcmc( as.numeric(chain )), prob = conf/100))[1]
  
  ucl <-as.numeric(HPDinterval( as.mcmc( as.numeric(chain ) ), prob = conf/100))[2]
  
  return(list(est=est,lcl=lcl,ucl=ucl))
  
  
}


frac.mixt.true <- function( mu , sigma , omega ,  c.oel) {
  
  
  res <-100*(1-pnorm((log((c.oel))-mu)/sigma))
  
  res <- res*(1-omega)
  
  return(res)
  
  
}

##### mixture and trunc  model estimation of the arithmetic mean


am.mixt <- function( mu.chain , sigma.chain , omega.chain ,  conf ) {
  
  
  chain <- exp(mu.chain + 0.5*sigma.chain^2)
  
  chain <- chain*(1-omega.chain)
  
  est <-median(chain)
  
  lcl <-as.numeric(HPDinterval( as.mcmc( as.numeric(chain )), prob = conf/100))[1]
  
  ucl <-as.numeric(HPDinterval( as.mcmc( as.numeric(chain ) ), prob = conf/100))[2]
  
  return(list(est=est,lcl=lcl,ucl=ucl))
  
  
}


am.mixt.true <- function( mu , sigma , omega ) {
  
  
  res <- exp(mu + 0.5*sigma^2)
  
  res <- res*(1-omega)
  
  return(res)
  
  
}




###### censored model 

perc.cens <- function( mu.chain , sigma.chain , target_perc = 0.5 , conf) {
  
  
  chain <-exp(mu.chain + qnorm(target_perc)*sigma.chain)
  
  est <-median(chain)
  
  lcl <-as.numeric(HPDinterval( as.mcmc( as.numeric(chain )), prob = conf/100))[1]
  
  ucl <-as.numeric(HPDinterval( as.mcmc( as.numeric(chain ) ), prob = conf/100))[2]
  
  return(list(est=est,lcl=lcl,ucl=ucl))
  
  
}


###### censored models  : estimation of the exceeedance 

frac.cens <- function( mu.chain , sigma.chain , conf , c.oel) {
  
  
  chain <-100*(1-pnorm((log((c.oel))-mu.chain)/sigma.chain))
  
  est <-median(chain)
  
  lcl <-as.numeric(HPDinterval( as.mcmc( as.numeric(chain )), prob = conf/100))[1]
  
  ucl <-as.numeric(HPDinterval( as.mcmc( as.numeric(chain ) ), prob = conf/100))[2]
  
  return(list(est=est,lcl=lcl,ucl=ucl))
  
  
}

###### censored models  : estimation of the arithmetic mean

am.cens <- function( mu.chain , sigma.chain , conf ) {
  
  
  chain <- exp( mu.chain + 0.5*sigma.chain^2)
  
  est <-median(chain)
  
  lcl <-as.numeric(HPDinterval( as.mcmc( as.numeric(chain )), prob = conf/100))[1]
  
  ucl <-as.numeric(HPDinterval( as.mcmc( as.numeric(chain ) ), prob = conf/100))[2]
  
  return(list(est=est,lcl=lcl,ucl=ucl))
  
  
}

####### wrapping function : GM/GSD/omega/P50/75/95+frac+am from chain from mixture and trunc models


fun.params.mixt <- function( stan.res , conf , c.oel) {
  
  
  return ( list( gm = gm( mu.chain = stan.res$mu , conf = conf) ,
                 
                 gsd = gsd( sigma.chain = stan.res$sigma , conf = conf) ,
                 
                 omega = omega.mixt( omega.chain = stan.res$omega , conf = conf) ,
    
                 med = perc.mixt( mu.chain = stan.res$mu , 
                                  sigma.chain = stan.res$sigma , 
                                  omega.chain = stan.res$omega,
                                  target_perc = 0.5 , conf = conf) ,
                 P75 = perc.mixt( mu.chain = stan.res$mu , 
                                  sigma.chain = stan.res$sigma , 
                                  omega.chain = stan.res$omega,
                                  target_perc = 0.75 , conf = conf),
                 
                 P95 = perc.mixt( mu.chain = stan.res$mu , 
                                  sigma.chain = stan.res$sigma , 
                                  omega.chain = stan.res$omega,
                                  target_perc = 0.95 , conf = conf),
                 
                 frac = frac.mixt( mu.chain = stan.res$mu , 
                                  sigma.chain = stan.res$sigma , 
                                  omega.chain = stan.res$omega ,
                                  conf = conf , c.oel = c.oel),
                 
                 am = am.mixt( mu.chain = stan.res$mu , 
                                 sigma.chain = stan.res$sigma , 
                                 omega.chain = stan.res$omega,
                                 conf = conf)
                 
                 ) )
  
  
  }

####### wrapping function : P50/75/95+frac+am from chain for censored model


fun.params.cens <- function( stan.res , conf , c.oel) {
  
  
  return ( list( gm = gm( mu.chain = stan.res$mu , conf = conf) ,
                 
                 gsd = gsd( sigma.chain = stan.res$sigma , conf = conf) ,
                 
                 omega = 0 ,
                 
                 med = perc.cens( mu.chain = stan.res$mu , 
                                  sigma.chain = stan.res$sigma ,
                                  target_perc = 0.5 , conf = conf) ,
                 P75 = perc.cens( mu.chain = stan.res$mu , 
                                  sigma.chain = stan.res$sigma ,
                                  target_perc = 0.75 , conf = conf),
                 
                 P95 = perc.cens( mu.chain = stan.res$mu , 
                                  sigma.chain = stan.res$sigma,
                                  target_perc = 0.95 , conf = conf),
                 
                 frac = frac.cens( mu.chain = stan.res$mu , 
                                   sigma.chain = stan.res$sigma ,
                                   conf = conf , c.oel = c.oel),
                 
                 am = am.cens( mu.chain = stan.res$mu , 
                                   sigma.chain = stan.res$sigma ,
                                   conf = conf)
                 
  ) )
  
  
}