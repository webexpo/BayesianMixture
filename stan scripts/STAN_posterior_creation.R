##########################################################################################
#
#    Function generating chains for the STAN MODELs for one dataset
#    
#    V0.6 July / 2026
#  
###########################################################################################


stan_bayes.mixt <- function( data.sample , n.iter , n.warmup , mymodel) {

####### PREPARING DATA FOR STAN

        # centered response
  
        centering <- median(data.sample$x)
  
        x <- data.sample$x /centering
  
        
        # logical vector identifying observed/unobserved
        
        is_observed <-as.integer(data.sample$x.is.observed)
        
        # List of input for the STAN models
        
        dataList = list( "x" = array(x,dim=length(x)) , 
                         "N" = data.sample$n,
                         "N0" = data.sample$n0,
                         "is_observed" = array(is_observed,dim=length(is_observed))
        )
        
        #Initial values : we leave them to STAN
        
      

###### GENERATING POSTERIOR SAMPLES

          # model initialization
          
          stanFit = sampling( object=mymodel , data=dataList, 
                            
                            chains=4 , 
                            
                            iter = n.iter  ,
                            
                            warmup = n.warmup , thin=1, show_messages=FALSE , control = list(adapt_delta = 0.8))
          

          
          ###### generating the MCMC chains
          
          #amples <-extract(stanFit,c("mu","sigma","omega"),permuted = TRUE, inc_warmup = FALSE)
          
          
########## results : centering point + chains + stanfit object
          
          return(list(stanFit = stanFit,
                 centering = centering))
            
            
}        
          



##################### CENSORED

stan_bayes.cens <- function( data.sample , n.iter , n.warmup , mymodel) {
  
  ####### PREPARING DATA FOR STAN
  
  # centered response
  
  centering <- median(data.sample$x)
  
  x <- data.sample$x /centering
  
  
  # logical vector identifying observed/unobserved
  
  is_observed <-as.integer(data.sample$x.is.observed)
  
  # List of input for the STAN models
  
  dataList = list( "x" = array(x,dim=length(x)) , 
                   "N" = data.sample$n,
                   "is_observed" = array(is_observed,dim=length(is_observed))
  )
  
  
  #initial values  - We leave them to STAN
  
  ###### GENERATING POSTERIOR SAMPLES
  
  # model initialization
  
  stanFit = sampling( object=mymodel , data=dataList, 
                      
                      chains=4 , 
                      
                      iter = n.iter  ,
                      
                      warmup = n.warmup , thin=1, show_messages=FALSE , control = list(adapt_delta = 0.8))
  
  
  
  ###### generating the MCMC chains
  
  #amples <-extract(stanFit,c("mu","sigma","omega"),permuted = TRUE, inc_warmup = FALSE)
  
  
  ########## results : centering point + chains + stanfit object
  
  return(list(stanFit = stanFit,
              centering = centering))
  
}        

