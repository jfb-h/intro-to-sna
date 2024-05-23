
# effective size

egonet.effsize <- function(egonet){
  n <- network.size(egonet)-1
  t <- network.edgecount(egonet)-n
  return(n-2*t/n)
}

network.effsize <- function(net){
  if (is.directed(net)) stop("Not implemented for directed graphs.")
  egonets <- lapply(ego.extract(net), FUN = network, directed = FALSE)
  return(sapply(nets, FUN = effsize))
}

net = flomarriage

net

# efficiency

egonet.efficiency <- function(net){

  efficiency <- function(egonet){
    n <- network.size(egonet)-1
    t <- network.edgecount(egonet)-degree(egonet, nodes = 1, cmode = "freeman",
                                          gmode = ifelse(is.directed(egonet),"digraph","graph"))
    efficiency <- 1-(2*t)/(n^2)
    return(efficiency)
  }
  nets <- lapply(ego.extract(net), FUN = network, directed = is.directed(net))
  return(sapply(nets, FUN = efficiency))
}



egonet.constraint <- function(net){
  
  rownorm <- function(x){
    rs <- rowSums(x)
    for(i in 1:nrow(x)){
      x[i,] <- x[i,]/sum(x[i,])
    }
    return(x)
  }
  
  constraint <- function(egonet){
    n <- nrow(egonet)-1
    
    egonet_rownorm <- rownorm(egonet)
    pij <- egonet_rownorm[1,-1]
    
    if(n == 0) return(NA)
    
    else if(n == 1) return(sum(pij))
    
    else{
      dc <- c()
      for(j in 1:n){
        piq <- pij[1:n != j]
        pqj <- egonet_rownorm[2:(n+1),2:(n+1)][1:n != j, j]
        dc[j] <- (pij[j] + c(piq %*% pqj))^2
      }
      return(sum(dc))
    }
  }
  egonets <- ego.extract(net)
  return(sapply(egonets, FUN = constraint))
}



egonet.hierarchy <- function(net){
  
  rownorm <- function(x){
    rs <- rowSums(x)
    for(i in 1:nrow(x)){
      x[i,] <- x[i,]/sum(x[i,])
    }
    return(x)
  }
  
  constraint <- function(egonet){
    n <- nrow(egonet)-1
    
    egonet_rownorm <- rownorm(egonet)
    pij <- egonet_rownorm[1,-1]
    
    if(n == 0) return(NA)
    else if(n == 1) return(sum(pij))
    else {
      dc <- c()
      for(j in 1:n){
        piq <- pij[1:n != j]
        pqj <- egonet_rownorm[2:(n+1),2:(n+1)][1:n != j, j]
        dc[j] <- (pij[j] + c(piq %*% pqj))^2
      }
      
      return(dc)
    }
  }
  
  hierarchy <- function(egonet){
    n <- nrow(egonet)-1
    if(n > 0){
      
      c <- constraint(egonet)
      conratio <- c/(sum(c)/n)
      num <- sum(conratio*log(conratio))
      den <- n*log(n)
      h <- ifelse(n == 1, 1, num/den)
      
      return(h)
    }
    else return(NA)
  }
  egonets <- ego.extract(net)
  return(sapply(egonets, FUN = hierarchy))
}
