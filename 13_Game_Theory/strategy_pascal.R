
strategyPascal1 <- function(ownpastactions, partnerpastactions){ # tit for tat
  
  n <- length(ownpastactions)
  i <- c(1:(n - 1))
  
  if (n == 0){return (T)}
  
  else {return (partnerpastactions[-i])}
    
}

stragetyPascal2 <- function(ownpastactions, partnerpastactions){ # tit for two tats
  
  n <- length(ownpastactions)
  i <- c(1:(n - 2))
  
  if (n == 0){return (T)}
  
  if (all(!partnerpastactions[-i])){
    return (F)
  
  else{return (T)}
}

boxp


if (c(T, F)){
  print("err")
}
