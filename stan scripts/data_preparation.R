#
#
#
#  DATA PREPARATION FOR THE BERNOUILLI LOGNORMAL MIXTURE MODEL WEB APPLICATION
#
# V0.5 July 12 / 2019
#


data.formatting.mixt <-function(data.in) {
  
  ##data.in is in the format of what is entered in expostats too1
  

  ######  Part A :  initial formating  
  
  ##finding the first digit
  
  
  Min <-regexpr('[0123456789]',data.in)
  
  if (substring(data.in,Min-1,Min-1)=='<') Min <-Min-1

  
  ##finding the last
  
  reva <-paste(rev(strsplit(data.in, split = "")[[1]]), collapse = "") 
  
  Max <-nchar(data.in)-regexpr('[0123456789]',reva)+1
  
  #cutting...
  
  data.in <-substring(data.in,Min,Max)
  
  #splitting into observations
  
  data.in <- strsplit(data.in,"\n")
  data.in <- unlist(data.in)
  data.in <-data.in[!is.na(data.in)]
  data.in <-data.in[!data.in==""]
  
  
  ######  Part B :  Information about censorship for further analysis 
  
  #preparing the result object
  
  result <-list()
  
  result$data <-data.in
  

  result$x.is.observed <-!grepl('<' , data.in, fixed = TRUE) 
  
  suppressWarnings(result$x <- as.numeric(data.in))
  
  result$x[ !result$x.is.observed ] <- as.numeric( substring( data.in[!result$x.is.observed] , 2) )
  
  result$n <- length(data.in)
  
  result$n0 <- length( data.in[ !result$x.is.observed ])

  result$LOQ <- result$x[!result$x.is.observed]
  

  return(result)
  
  
}

######################## data formatting for the EXCEL format

data.formatting.file <-function(myfile) {
  
  ##data.in is in the format of what is entered in expostats too1
  
  
  data.in <-unlist(myfile[,1])
  
  #preparing the result object
  
  result <-list()
  
  result$data <-data.in
  
  result$x.is.observed <-!grepl('<' , data.in, fixed = TRUE) 
  
  suppressWarnings(result$x <- as.numeric(data.in))
  
  result$x[ !result$x.is.observed ] <- as.numeric( substring( data.in[!result$x.is.observed] , 2) )
  
  result$n <- length(data.in)
  
  result$n0 <- length( data.in[ !result$x.is.observed ])
  
  result$LOQ <- result$x[!result$x.is.observed]
  
  
  return(result)
  
  
}

