######################################
#
#  FINAL STAN MODELS
#  
#  STAN model statements
#
#  V0.5 September 2019
#  V0.6 July 2026 : new STAN syntax updated for array declarations
#   NOTE : The models need to be recompiled on each system
######################################

library(rstan)

####################################  MIXTURE MODEL

model.stan.mixt.text  = "
          
          data{
          int<lower=0> N;
          int<lower=0> N0;
          array[N] real<lower=0> x;
          array[N] int<lower=0,upper=1> is_observed;

          }
          
          parameters{
          real mu;
          real logsigma;
          real<lower=0,upper=1> Pmax;
          real<lower=0,upper=Pmax> omega;
          }
          
          transformed parameters {

          real <lower=0> sigma;
        
          sigma = exp(logsigma);
          
          }
          
          model{
          
          //priors
          
          mu~normal(0,3.16);
          
          logsigma~normal(-0.1744,0.6259421);
          
          Pmax~beta( N0 , ( N - N0 ) );
          
          omega~uniform( 0 , Pmax );
          
          
          //likelihood
          
          for (n in 1:N) {
          
          if (is_observed[n]==0)   // For censored values
          
          target += log_sum_exp( bernoulli_lpmf( 1 | omega ) ,  bernoulli_lpmf(0 | omega ) + lognormal_lcdf( x[n] | mu, sigma) ) ;
          
          else     // for observed values, classic lognormal probability density
          
          target += bernoulli_lpmf( 0 | omega ) + lognormal_lpdf( x[n] | mu, sigma);
          
          }
          
          }
          
          //  extracting the likelihood
          generated quantities {
            
            vector[N] log_lik;
            
            for (n in 1:N) {
            
          if (is_observed[n]==0)   // For censored values
          
         log_lik[n] = log_sum_exp( bernoulli_lpmf( 1 | omega ) ,  bernoulli_lpmf(0 | omega ) + lognormal_lcdf( x[n] | mu, sigma) ) ;
          
          else     // for observed values, classic lognormal probability density
          
          log_lik[n] = bernoulli_lpmf( 0 | omega ) + lognormal_lpdf( x[n] | mu, sigma);
       
            }
          
          }
          
          " # close quote for modelString


stanmodel.mixt <- stan_model( model_code = model.stan.mixt.text  )


####################################  CENSORED MODEL



stanmodel.cens.text = "
          
          data{
          int<lower=0> N;
          array[N] real<lower=0> x;
          array[N] int<lower=0,upper=1> is_observed;
          }
          
          parameters{
          real mu;
          real logsigma;
          }
          
          transformed parameters {
          real <lower=0> sigma;
          sigma = exp(logsigma);
          }
          
          model{
          
          //priors
          
          mu~normal(0,3.16);
          
          logsigma~normal(-0.1744,0.6259421);
          
          
          //likelihood
          
          for (n in 1:N) {
          
          if (is_observed[n]==0)   // For censored values
          
          target += lognormal_lcdf( x[n] | mu, sigma);  
          
          else     // for observed values
          
          target += lognormal_lpdf( x[n] | mu, sigma);
          
          }
          
          }
          
        //  extracting the likelihood
          generated quantities {
            
            vector[N] log_lik;
            
          for (n in 1:N) {
          
          if (is_observed[n]==0)   // For censored values
          
         log_lik[n] = lognormal_lcdf( x[n] | mu, sigma) ; 
          
          else     // for observed values
          
          log_lik[n] = lognormal_lpdf( x[n] | mu, sigma);
          
          }
          
          }
          " # close quote for modelString


stanmodel.cens <- stan_model( model_code = stanmodel.cens.text )





