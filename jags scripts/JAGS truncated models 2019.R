############################
#
#  truncated model in JAGS 
#
#
################################


###########

model.jags.A <-paste("

model { 
  
  ### Priors
  
      mu ~ dnorm(0, 0.01) #vague on log-mean
      
      sigma <-exp(log.sigma)
          
      log.sigma ~ dnorm(-0.1744,2.5523)
      
  ## prior for proportion of true zeroes    

      Pmax ~ dbeta( n.unobs , (n.total-n.unobs) ) 
    
      P.true0 ~ dunif( 0 , Pmax )
  
  #likelihood
  #about detected above Xd that tell us about mu and sigma
  #independent of anything that goes on below Xd
  #uses truncated normal
      
      
      for (i in 1:n.obs) 
              
                     {
              
              X[i] ~ dlnorm(mu, tau)T(Xd, ) 
              
              } 
      #priors
      mu ~ dnorm(0, 0.01) #vague on log-mean
      tau~dgamma(a,b) # precision: fix a and b to match Sym-Kromhout on log-st.dev
      a<-0.1
      b<-0.1
      sigma <- 1/sqrt(tau) 
      Xd<- 0.4386971 #can be changed to vary and depend on duration of sampling T.
  
  #calculation MLE by Taylor
      epsilon<-(log(Xd)-mu)/sigma
      omega<-(n0-n*phi(epsilon)) / (n-n*phi(epsilon))
      
}

 
                            ",sep='')