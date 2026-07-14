############################
#
#  bernoulli lognormal mixture model in JAGS using the zeroes trick
#
#
################################

####### A. Probit PRIOR
            
            model.jags.A <-paste("
                            
                            model {
                            
                            
                            ##priors for the lognormal part (expostats default priors)
                            
                                mu  ~dunif(-100,100)
                                
                                sigma <-exp(log.sigma)
                                
                                log.sigma ~ dnorm(-0.1744,2.5523)
                                
                            ## prior for proportion of true zeroes    

                                P.true0 <-ilogit(alpha)
                                
                                alpha ~ dnorm( 0 , 0.1 )
                                
                            
                            ###likelihood
                            
                            for (i in 1:n.total) {
                            
                                # the zeroes tricks : allows to write the log-likelihood in JAGS (ll[i])
                                
                                zeros[i] ~ dpois(-ll[i] + C)
                                
                            
                                # log-likelihood for the lognormal distribution
                                
                                
                                ln1[i] <- -(log(y[i]) + log(sigma) - log(sqrt(2 * pi.constant)))
                                
                                ln2[i] <- -0.5 * pow((log(y[i]) - mu),2)/(sigma * sigma)
                                
                                LN[i] <- ln1[i] + ln2[i]
                                
                            
                                #mixture likelihood
                                
                                ll[i] <- ifelse(is.observed[i],
                                            
                                           log(1-P.true0) + LN[i] ,   #### for observed
    
                                            log(P.true0 + (1-P.true0)*phi( (log(y[i])-mu)/sigma ) ) ###for not observed (zero or censored)
    
                                    )
                            
                            
                            }
                            

                            }
                            
                            ",sep='')
            
            
            
####### B. Uniform prior
            
            model.jags.B <-paste("
                            
                            model {
                            
                            
                            ##priors for the lognormal part (expostats default priors)
                            
                                mu  ~dunif(-100,100)
                                
                                sigma <-exp(log.sigma)
                                
                                log.sigma ~ dnorm(-0.1744,2.5523)
                                
                            ## prior for proportion of true zeroes    

                                P.true0 ~ dunif(0,1)
                                
                            
                            ###likelihood
                            
                            for (i in 1:n.total) {
                            
                                # the zeroes tricks : allows to write the log-likelihood in JAGS (ll[i])
                                
                                zeros[i] ~ dpois(-ll[i] + C)
                                
                            
                                # log-likelihood for the lognormal distribution
                                
                                
                                ln1[i] <- -(log(y[i]) + log(sigma) - log(sqrt(2 * pi.constant)))
                                
                                ln2[i] <- -0.5 * pow((log(y[i]) - mu),2)/(sigma * sigma)
                                
                                LN[i] <- ln1[i] + ln2[i]
                                
                            
                                #mixture likelihood
                                
                                ll[i] <- ifelse(is.observed[i],
                                            
                                           log(1-P.true0) + LN[i] ,   #### for observed
    
                                            log(P.true0 + (1-P.true0)*phi( (log(y[i])-mu)/sigma ) ) ###for not observed (zero or censored)
    
                                    )
                            
                            
                            }
                            

                            }
                            
                            ",sep='')
            
####### C. Uniform prior - Igors' proposal
            
            model.jags.C <-paste("
                            
                            model {
                            
                            
                            ##priors for the lognormal part (expostats default priors)
                            
                                mu  ~dunif(-100,100)
                                
                                sigma <-exp(log.sigma)
                                
                                log.sigma ~ dnorm(-0.1744,2.5523)
                                
                            ## prior for proportion of true zeroes 
                            
                                Pmax ~ dbeta( n.unobs , (n.total-n.unobs) ) 

                                P.true0 ~ dunif( 0 , Pmax )
                                
                                
                                
                            
                            ###likelihood
                            
                            for (i in 1:n.total) {
                            
                                # the zeroes tricks : allows to write the log-likelihood in JAGS (ll[i])
                                
                                zeros[i] ~ dpois(-ll[i] + C)
                                
                            
                                # log-likelihood for the lognormal distribution
                                
                                
                                ln1[i] <- -(log(y[i]) + log(sigma) - log(sqrt(2 * pi.constant)))
                                
                                ln2[i] <- -0.5 * pow((log(y[i]) - mu),2)/(sigma * sigma)
                                
                                LN[i] <- ln1[i] + ln2[i]
                                
                            
                                #mixture likelihood
                                
                                ll[i] <- ifelse(is.observed[i],
                                            
                                           log(1-P.true0) + LN[i] ,   #### for observed
    
                                            log(P.true0 + (1-P.true0)*phi( (log(y[i])-mu)/sigma ) ) ###for not observed (zero or censored)
    
                                    )
                            
                            
                            }
                            

                            }
                            
                            ",sep='')
            
            