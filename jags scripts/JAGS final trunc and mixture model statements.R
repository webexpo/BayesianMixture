#############################################################
#
#   statements for the mixture and truncated models
#
#   These statement are reference ; they shouls be used in preference to all other versions
#
###########################################################


######### MIXTURE MODEL

model.jags.mixture <-paste("
                                 
                                 model {
                                 
                                 
                                 ####priors for the lognormal part 
                                 
                                 mu ~ dnorm( 0 , 0.1)              #vague on log-mean
                                 
                                 sigma <-exp(log.sigma)            # weakly informative: based on historical variability data 
                                                                   # described in Lavoue et al., 2018, Ann. Work. Expo. Health 63(3):267–279
                                 log.sigma ~ dnorm(-0.1744,2.5523)
                                 
                                 ## prior for proportion of true zeroes (omega)
                                 
                                 Pmax ~ dbeta( n0 , ( n - n0 ) )    # prior on maximum value of omega (omega_max in text)
                                 
                                 omega ~ dunif( 0 , Pmax )          # prior on omega  
                                 
                                 
                                 ####likelihood function 
                                 
                                 for (i in 1:n) {
                                 
                                 # the zeroes tricks : allows to write the log-likelihood in JAGS (ll[i])
                                 
                                 zeros[i] ~ dpois(-ll[i] + C)
                                 
                                 
                                 # log-likelihood for the lognormal distribution
                                 
                                 
                                 LN[i] <- -(log(x[i]) + log(sigma) - log(sqrt(2 * pi.constant))) -0.5 * pow((log(x[i]) - mu),2)/(sigma * sigma)
                                 
                                 
                                 #mixture likelihood : 2 different expressions depending on observed or not
                                 
                                 ll[i] <- ifelse(is.observed[i],
                                 
                                 # for observed values : likelihood is (1-omega)*lognormal likelihood
                                 
                                         log(1-omega) + LN[i] ,   
                                 
                                 # for unobserved values : likelihood is part point mass at 0 (with probability omega),
                                 # part censored lognormal (using the cumulative density), with probability (1-omega)
                                 
                                        log(omega + (1-omega)*phi( (log(x[i])-mu)/sigma ) ) 
                                 
                                 )
                                 
                                 
                                 }
                                 
                                 
                                 }
                                 
                                 ",sep='')  




####### TRUNCATED MODEL

model.jags.truncated <- paste("        
                                    
                                    model { 
                                    
                                    ### priors          
                                    
                                    mu ~ dnorm( 0 , 0.1)              #vague on log-mean
                                    
                                    tau <- 1/(sigma^2)
                                    
                                    sigma <-exp(log.sigma)            # weakly informative based on historical variability data 
                                                                      # described in Lavoue et al., 2018, Ann. Work. Expo. Health 63(3):267–279
                                    log.sigma ~ dnorm(-0.1744,2.5523)
                                    
                                    
                                    
                                    #likelihood
                                    
                                    #only detected values above Xd tell us about mu and sigma
                                    #independent of anything that goes on below Xd
                                    #uses truncated lognormal
                                    
                                    
                                    for (i in 1:length(X1) ) {X1[i] ~ dlnorm(mu, tau)T(Xd, ) } 
                                    
                                    #calculation of epsilon and omega in Taylor et al. from the estimated truncated distribution 
                                    
                                    epsilon <- (log(Xd)-mu)/sigma
                                    
                                    omega <- (n0-n*phi(epsilon)) / (n-n*phi(epsilon))
                                    
                                    }
                                    
                                    ",sep='') 