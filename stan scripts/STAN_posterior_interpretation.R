#####################################################
#
#   results functions for the MIXT / CENS models
#
#    V0.5 Sept / 2019
#    V0.6 July / 2026 partial implemeentation of web app scripts (no graphs / posterior predictive checks)
#
#######################################################


#################### chains

posterior.chains <- function( mixture.analysis , censored.analysis  ) {
  
  mix <- extract(mixture.analysis$stanFit,c("mu","sigma","omega"),permuted = TRUE, inc_warmup = FALSE)
  
  cens <- extract(censored.analysis$stanFit,c("mu","sigma"),permuted = TRUE, inc_warmup = FALSE)
  
  
  return(list(
  
  gm.mix = as.numeric(exp( mix$mu ) * mixture.analysis$centering) ,
  gm.cens = as.numeric(exp( cens$mu ) * censored.analysis$centering) ,
  
  gsd.mix = as.numeric(exp( mix$sigma )) ,
  gsd.cens = as.numeric(exp( cens$sigma )) ,
  
  omega.mix = mix$omega ))

  }
  
#################### table of point estimates and HPDIs


estimates <- function( mixture.analysis , censored.analysis , conf ) {
  
  
  
  results <- data.frame( model = c("Mixture" , "Censored") ,
                         GM = character(2) ,
                         GSD = character(2) ,
                         omega = c( "" , "-" ) , stringsAsFactors = FALSE)
  
  mix <- extract(mixture.analysis$stanFit,c("mu","sigma","omega"),permuted = TRUE, inc_warmup = FALSE)
  
  cens <- extract(censored.analysis$stanFit,c("mu","sigma"),permuted = TRUE, inc_warmup = FALSE)
  
  gm.chain.mix <- as.numeric( exp( mix$mu ) * mixture.analysis$centering )
  GM.mix.est <- as.numeric(median( gm.chain.mix ))
  GM.mix.int <- as.numeric(HPDinterval( as.mcmc( gm.chain.mix ), prob = conf/100))
  
  gm.chain.cens <- as.numeric( exp( cens$mu ) * censored.analysis$centering )
  GM.cens.est <- as.numeric(median( gm.chain.cens ))
  GM.cens.int <- as.numeric(HPDinterval( as.mcmc( gm.chain.cens ), prob = conf/100))
  
  results$GM[1] <- paste( signif( GM.mix.est , 2 ),
                          " [ ",
                          signif( GM.mix.int[1] ,2 ),
                          " - ",
                          signif( GM.mix.int[2] ,2 ),
                          " ]" , sep="")
  
  
  results$GM[2] <- paste( signif( GM.cens.est , 2 ),
                          " [ ",
                          signif( GM.cens.int[1] ,2 ),
                          " - ",
                          signif( GM.cens.int[2] ,2 ),
                          " ]" , sep="")
  
  
  gsd.chain.mix <- as.numeric( exp( mix$sigma ) )
  gsd.mix.est <- as.numeric(median( gsd.chain.mix ))
  gsd.mix.int <- as.numeric(HPDinterval( as.mcmc( gsd.chain.mix ), prob = conf/100))
  
  gsd.chain.cens <- as.numeric( exp( cens$sigma ) )
  gsd.cens.est <- as.numeric(median( gsd.chain.cens ))
  gsd.cens.int <- as.numeric(HPDinterval( as.mcmc( gsd.chain.cens ), prob = conf/100))
  
  results$GSD[1] <- paste( signif( gsd.mix.est , 2 ),
                           " [ ",
                           signif( gsd.mix.int[1] ,2 ),
                           " - ",
                           signif( gsd.mix.int[2] ,2 ),
                           " ]" , sep="")
  
  
  results$GSD[2] <- paste( signif( gsd.cens.est , 2 ),
                           " [ ",
                           signif( gsd.cens.int[1] ,2 ),
                           " - ",
                           signif( gsd.cens.int[2] ,2 ),
                           " ]" , sep="")
  
  
  omega.chain.mix <- as.numeric( mix$omega )
  omega.mix.est <- as.numeric(median( omega.chain.mix ))
  omega.mix.int <- as.numeric(HPDinterval( as.mcmc( omega.chain.mix ), prob = conf/100))
  
  results$omega[1] <- paste( signif( omega.mix.est , 2 ),
                             " [ ",
                             signif( omega.mix.int[1] ,2 ),
                             " - ",
                             signif( omega.mix.int[2] ,2 ),
                             " ]" , sep="")
  
  return(results)
  
}


#################### goodness of fit


goodness.of.fit <- function( mixture.analysis , censored.analysis ) {
  
  
  mix <- mixture.analysis$stanFit
  
  cens <- censored.analysis$stanFit
  
  #### recuparating the log-likelhoods at each point / iteration, and calculating the WAICs
  
  log_lik.mixt <- extract_log_lik( stanfit = mix, merge_chains = TRUE)
  
  waic.mix <-waic(log_lik.mixt)$estimates[3,1]
  
  log_lik.cens <- extract_log_lik( stanfit = cens, merge_chains = TRUE)
  
  waic.cens <-waic(log_lik.cens)$estimates[3,1]
  
  
  ### calculating weights and the odds of a better mixture model
  
  vec.waic <- c( waic.mix ,waic.cens)
  
  denominator <- sum( exp( -0.5 * ( vec.waic - min( vec.waic )) ) ) 
  
  weight.mix <- exp( - 0.5 * ( waic.mix - min( vec.waic ) ) ) / denominator
  
  weight.cens <- exp( - 0.5 * ( waic.cens - min( vec.waic ) ) ) / denominator
  
  Odds <- weight.mix/weight.cens
  
  
  ### final results
  
  results <-list( weight.mix = signif(weight.mix,2) ,
                  weight.cens = signif(weight.cens,2) ,
                  Odds.mixovercens = signif(Odds,2),
                  waic.mix =  waic.mix,
                  waic.cens =  waic.cens)
  
  return(results)
  
}


#################### table of point other metrics and HPDIs


othermetrics <- function( mixture.analysis , censored.analysis , conf , oel , target_perc ) {
  
  
  
  results <- data.frame( Model = c("Mixture" , "Censored") ,
                         Percentile = character(2) ,
                         Exceedance = character(2) ,
                         AM = character(2) , stringsAsFactors = FALSE)
  
  mix <- extract(mixture.analysis$stanFit,c("mu","sigma","omega"),permuted = TRUE, inc_warmup = FALSE)
  
  cens <- extract(censored.analysis$stanFit,c("mu","sigma"),permuted = TRUE, inc_warmup = FALSE)
  
  # percentile
  
  perc.chain.mix <- perc.mixt.chain( mu.chain = mix$mu + log(mixture.analysis$centering), sigma.chain = mix$sigma , omega.chain = mix$omega , target_perc = target_perc , conf = conf)
  perc.mix.est <- as.numeric(median( perc.chain.mix ))
  perc.mix.int <- as.numeric(HPDinterval( as.mcmc( perc.chain.mix ), prob = conf/100))
  
  perc.chain.cens <- perc.cens.chain( mu.chain = cens$mu + log(censored.analysis$centering), sigma.chain = cens$sigma , target_perc = target_perc , conf = conf)
  perc.cens.est <- as.numeric(median( perc.chain.cens ))
  perc.cens.int <- as.numeric( HPDinterval( as.mcmc( perc.chain.cens ), prob = conf/100))
  
  results$Percentile[1] <- paste( signif( perc.mix.est , 2 ),
                          " [ ",
                          signif( perc.mix.int[1] ,2 ),
                          " - ",
                          signif( perc.mix.int[2] ,2 ),
                          " ]" , sep="")
  
  
  results$Percentile[2] <- paste( signif( perc.cens.est , 2 ),
                          " [ ",
                          signif( perc.cens.int[1] ,2 ),
                          " - ",
                          signif( perc.cens.int[2] ,2 ),
                          " ]" , sep="")
  
  # Exceedance
  
  
  frac.chain.mix <- frac.mixt.chain( mu.chain = mix$mu + log(mixture.analysis$centering), sigma.chain = mix$sigma , omega.chain = mix$omega , conf = conf , c.oel = oel)
  frac.mix.est <- as.numeric(median( frac.chain.mix ))
  frac.mix.int <- as.numeric(HPDinterval( as.mcmc( frac.chain.mix ), prob = conf/100))
  
  frac.chain.cens <- frac.cens.chain( mu.chain = cens$mu + log(censored.analysis$centering), sigma.chain = cens$sigma , conf = conf , c.oel = oel)
  frac.cens.est <- as.numeric(median( frac.chain.cens ))
  frac.cens.int <- as.numeric(HPDinterval( as.mcmc( frac.chain.cens ), prob = conf/100))
  
  results$Exceedance[1] <- paste( signif( frac.mix.est , 2 ),
                                  " [ ",
                                  signif( frac.mix.int[1] ,2 ),
                                  " - ",
                                  signif( frac.mix.int[2] ,2 ),
                                  " ]" , sep="")
  
  
  results$Exceedance[2] <- paste( signif( frac.cens.est , 2 ),
                                  " [ ",
                                  signif( frac.cens.int[1] ,2 ),
                                  " - ",
                                  signif( frac.cens.int[2] ,2 ),
                                  " ]" , sep="")
  
  
  #### arithmetic mean
  
  
  am.chain.mix <- am.mixt.chain( mu.chain = mix$mu + log(mixture.analysis$centering), sigma.chain = mix$sigma , omega.chain = mix$omega , conf = conf )
  am.mix.est <- as.numeric(median( am.chain.mix ))
  am.mix.int <- as.numeric(HPDinterval( as.mcmc( am.chain.mix ), prob = conf/100))
  
  am.chain.cens <- am.cens.chain( mu.chain = cens$mu + log(censored.analysis$centering), sigma.chain = cens$sigma , conf = conf )
  am.cens.est <- as.numeric(median( am.chain.cens ))
  am.cens.int <- as.numeric(HPDinterval( as.mcmc( am.chain.cens ), prob = conf/100))
  
  results$AM[1] <- paste( signif( am.mix.est , 2 ),
                                  " [ ",
                                  signif( am.mix.int[1] ,2 ),
                                  " - ",
                                  signif( am.mix.int[2] ,2 ),
                                  " ]" , sep="")
  
  
  results$AM[2] <- paste( signif( am.cens.est , 2 ),
                                  " [ ",
                                  signif( am.cens.int[1] ,2 ),
                                  " - ",
                                  signif( am.cens.int[2] ,2 ),
                                  " ]" , sep="")
  
  return(results)
  
}



##############  Minor secondary functions


perc.mixt.chain <- function( mu.chain , sigma.chain , omega.chain , target_perc = 0.5 , conf) {
  
  
  quant <- ( target_perc - omega.chain) / ( 1 - omega.chain)
  
  chain <- numeric(length(mu.chain))
  
  criterion <- omega.chain >= target_perc
  
  chain[ criterion ] <- 0
  
  chain[ !criterion ] <- exp( qnorm( quant[ !criterion  ] , mu.chain[ !criterion  ] , sigma.chain[ !criterion  ] ) )
  
  return(as.numeric(chain))
  
  
}



perc.cens.chain <- function( mu.chain , sigma.chain , target_perc = 0.5 , conf) {
  
  
  chain <-exp(mu.chain + qnorm(target_perc)*sigma.chain)
  
  return(as.numeric(chain))
  
  
}


frac.mixt.chain <- function( mu.chain , sigma.chain , omega.chain ,  conf , c.oel) {
  
  
  chain <-100*(1-pnorm((log((c.oel))-mu.chain)/sigma.chain))
  
  chain <- chain*(1-omega.chain)
  
   return(as.numeric(chain))
  
  
}

frac.cens.chain <- function( mu.chain , sigma.chain , conf , c.oel) {
  
  
  chain <-100*(1-pnorm((log((c.oel))-mu.chain)/sigma.chain))
  
  return(as.numeric(chain))
  
  
}

am.mixt.chain <- function( mu.chain , sigma.chain , omega.chain ,  conf ) {
  
  
  chain <- exp(mu.chain + 0.5*sigma.chain^2)
  
  chain <- chain*(1-omega.chain)

  return(as.numeric(chain))
  
  
}


am.cens.chain <- function( mu.chain , sigma.chain , conf ) {
  
  
  chain <- exp( mu.chain + 0.5*sigma.chain^2)
  
  return(as.numeric(chain))
  
  
}
