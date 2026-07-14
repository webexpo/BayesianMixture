############ COMPARISON OF waic FOR THE 2 MODELS



######################################################################################


fun.compare.waic <-function(data0 , n.iter = 5000 , n.warmup = 2500 , n.chain = 1) {
  
  data.sample <-data0
  
  #### mixture model 
  
  centering <- median(data.sample$x)
  
  x <-data.sample$x / centering
  
  is_observed <-as.integer(data.sample$x.is.observed)
  
  dataList.mixt = list( "x" = array(x,dim=length(x)) , 
                   "N" = data.sample$n,
                   "N0" = data.sample$n0,
                   "is_observed" = array(is_observed,dim=length(is_observed))
  )
  
  
  
  stanFit.mixt = sampling( object=stanmodel.mixt , data=dataList.mixt , 
                           chains=n.chain , iter=n.iter , warmup=n.warmup , thin=1, show_messages=FALSE )
  
  samples.mixt <-extract(stanFit.mixt,c("mu","sigma","omega"),permuted = TRUE, inc_warmup = FALSE)
  
  
  # Extract pointwise log-likelihood and compute LOO
  log_lik.mixt <- extract_log_lik( stanfit = stanFit.mixt, merge_chains = FALSE)
  
  # as of loo v2.0.0 we can optionally provide relative effective sample sizes
  # when calling loo, which allows for better estimates of the PSIS effective
  # sample sizes and Monte Carlo error
  #r_eff.mixt <- relative_eff(exp(log_lik.mixt)) 
  
  #loo.mixt <- loo(log_lik.mixt, r_eff = r_eff.mixt, cores = 2)
  
  waic.mixt <-waic(log_lik.mixt)
  
  ### cleaning
  
  rm(stanFit.mixt)
  rm(log_lik.mixt)
  
  
  gc()
  
  ## censored model
  
  centering <- median(data.sample$x)
  
  x <-data.sample$x / centering
  
  is_observed <-as.integer(data.sample$x.is.observed)
  
  dataList.cens = list( "x" = array(x,dim=length(x)) , 
                        "N" = data.sample$n,
                        "is_observed" = array(is_observed,dim=length(is_observed))
  )
  
  
  stanFit.cens = sampling( object=stanmodel.cens , data=dataList.cens , 
                           chains=n.chain , iter=n.iter , warmup=n.warmup ,thin=1, show_messages=FALSE )
  
  
  samples.cens <-extract(stanFit.cens,c("mu","sigma"),permuted = TRUE, inc_warmup = FALSE)
  

  
  # Extract pointwise log-likelihood and compute LOO
  log_lik.cens <- extract_log_lik( stanfit = stanFit.cens, merge_chains = FALSE)
  
  # as of loo v2.0.0 we can optionally provide relative effective sample sizes
  # when calling loo, which allows for better estimates of the PSIS effective
  # sample sizes and Monte Carlo error
  #r_eff.cens <- relative_eff(exp(log_lik.cens)) 
  
  #loo.cens <- loo(log_lik.cens, r_eff = r_eff.cens, cores = 2)
  
  waic.cens <-waic(log_lik.cens)
  
  ### cleaning
  
  rm(stanFit.cens)
  rm(log_lik.cens)

  gc()
  
  ####
  
  
  ### comparison of the 2 models      
  
  
  #result.loo <- loo_compare(loo.mixt,loo.cens)
  
  #result.waic <- compare(waic.mixt,waic.cens)
  

  # model weights
  
  waic.vec <- c( waic.mixt$estimates[3,1] , waic.cens$estimates[3,1] )
  
  denominator <- sum( exp( -0.5 * (   waic.vec - min(waic.vec) ) ) )
  
  weight.mixt <- exp( - 0.5 * ( waic.vec[1] - min(waic.vec) ) ) / denominator
  
  weight.cens <- exp( - 0.5 * ( waic.vec[2] - min(waic.vec) ) ) / denominator
  
  
  OR.mixtoncens <-weight.mixt / weight.cens  
  
  ## comparison of estimates of MG and GSD (89% CI)
  
  res.table <-data.frame( model = c("mixture","censored"),
                          gm=character(2),
                          gsd=character(2), stringsAsFactors = FALSE)
  
      #mixture model
      
 
      mixt.chain.gm <-as.numeric(exp(samples.mixt$mu)*centering)
      
      mixt.hpdi.gm <- HPDinterval(  mcmc(mixt.chain.gm) , 0.89 )
      
      mixt.median.gm <-median(mixt.chain.gm)
      
      res.gm.mixt <- paste( signif( mixt.median.gm , 2 ) , " [" , signif( mixt.hpdi.gm[1] , 2) , " - " , signif( mixt.hpdi.gm[2] , 2) , "]" , sep="") 
      
      
      mixt.chain.gsd <-as.numeric(exp(samples.mixt$sigma))
      
      mixt.hpdi.gsd <- HPDinterval(  mcmc(mixt.chain.gsd) , 0.89 )
      
      mixt.median.gsd <-median(mixt.chain.gsd)
      
      res.gsd.mixt <- paste( signif( mixt.median.gsd , 2 ) , " [" , signif( mixt.hpdi.gsd[1] , 2) , " - " , signif( mixt.hpdi.gsd[2] , 2) , "]" , sep="") 
      
      #censored model
      

      cens.chain.gm <-as.numeric(exp(samples.cens$mu)*centering)
      
      cens.hpdi.gm <- HPDinterval(  mcmc(cens.chain.gm) , 0.89 )
      
      cens.median.gm <-median(cens.chain.gm)
      
      res.gm.cens <- paste( signif( cens.median.gm , 2 ) , " [" , signif( cens.hpdi.gm[1] , 2) , " - " , signif( cens.hpdi.gm[2] , 2) , "]" , sep="") 
      
      
      cens.chain.gsd <-as.numeric(exp(samples.cens$sigma))
      
      cens.hpdi.gsd <- HPDinterval(  mcmc(cens.chain.gsd) , 0.89 )
      
      cens.median.gsd <-median(cens.chain.gsd)
      
      res.gsd.cens <- paste( signif( mixt.median.gsd , 2 ) , " [" , signif( mixt.hpdi.gsd[1] , 2) , " - " , signif( mixt.hpdi.gsd[2] , 2) , "]" , sep="") 
      
      res.table$gm[1] <- res.gm.mixt
      res.table$gsd[1] <- res.gsd.mixt
      res.table$gm[2] <- res.gm.cens
      res.table$gsd[2] <- res.gsd.cens
      
      
      
  ## omega
  
      mixt.chain.omega <-as.numeric(samples.mixt$omega)
      
      mixt.hpdi.omega <- HPDinterval(  mcmc(mixt.chain.omega) , 0.89 )
      
      mixt.median.omega <-median(mixt.chain.omega)
      
      
     res.omega <- paste( signif( mixt.median.omega , 2 ) , " [" , signif( mixt.hpdi.omega[1] , 2) , " - " , signif( mixt.hpdi.omega[2] , 2) , "]" , sep="") 
  
     res.n0.n <- data.sample$n0/data.sample$n
     
     
  return(list(
    
    weight.mixt = weight.mixt,
    weight.cens = weight.cens,
    OR.mixtoncens = OR.mixtoncens,
    estimates = res.table,
    omega = res.omega,
    n0n = res.n0.n
    
      ))
}